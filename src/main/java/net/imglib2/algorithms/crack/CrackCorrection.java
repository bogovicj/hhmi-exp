package net.imglib2.algorithms.crack;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.ejml.data.DenseMatrix64F;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import edu.jhu.ece.iacl.algorithms.MGDM.MgdmDecomposition;
import edu.jhu.ece.iacl.utility.ArrayUtil;
import ij.IJ;
import ij.ImagePlus;
import net.imglib2.*;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.img.imageplus.*;
import net.imglib2.algorithm.edge.*;
import net.imglib2.algorithms.geometry.RodriguesRotation;
import net.imglib2.algorithms.moments.ImgMoment;
import net.imglib2.img.*;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.io.ImgIOException;
import net.imglib2.io.ImgOpener;
import net.imglib2.iterator.IntervalIterator;
import net.imglib2.ops.operation.randomaccessibleinterval.unary.DistanceMap;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.InverseRealTransform;
import net.imglib2.realtransform.InvertibleRealTransform;
import net.imglib2.realtransform.RealTransformRandomAccessible;
import net.imglib2.realtransform.RealTransformRandomAccessible.RealTransformRandomAccess;
import net.imglib2.realtransform.RealTransformRealRandomAccessible.RealTransformRealRandomAccess;
import net.imglib2.realtransform.RealViews;
import net.imglib2.type.numeric.*;
import net.imglib2.type.numeric.integer.*;
import net.imglib2.type.numeric.real.*;
import net.imglib2.view.Views;
import net.imglib2.view.composite.CompositeIntervalView;
import net.imglib2.view.composite.GenericComposite;
import net.imglib2.util.ImgUtil;
import net.imglib2.util.Util;

/**
 *
 * Class performing correction of high-pressure freezing cracks
 * in electron microscopy images.
 * 
 * Notes on implementation:
 * {@link Edgel}s are computed such that their position is <i>inside</i> the 
 * crack mask.  Local depth can be computed either by directly polling the distance
 * map at patch coordinates or by tracing a patch in the edgel normal direction until
 * a point outside the crack is reached.
 *
 * @author John Bogovic
 *
 * @param <T> image type
 * @param <B> label type
 */
public class CrackCorrection<T extends RealType<T>, B extends RealType<B>> {

	Img<T> img;
	Img<B> mask;

	Img<ByteType> edgelPatchMasks;
	byte b = 1;

	RealRandomAccess<T> interpRa;

	//IntegralImgDouble<B> maskIntImg;
	ArrayList<Edgel> edgels;
	
	Img<FloatType> depthPatch;
	Img<T> imgPatch;

	Img<B> imgChunks;

	private int maxNumNormalSteps = 30;
	
	protected static Logger logger = LogManager.getLogger(CrackCorrection.class
			.getName());

	public CrackCorrection() {
	}

	public CrackCorrection(Img<T> img, Img<B> mask) {
		this.img = img;
		this.mask = mask;

		NLinearInterpolatorFactory<T> interpFactory = new NLinearInterpolatorFactory<T>();

		RealRandomAccessible<T> interpolant = Views.interpolate(
				Views.extendMirrorSingle(img), interpFactory);

		interpRa = interpolant.realRandomAccess();

		ArrayImgFactory<ByteType> bfactory = new ArrayImgFactory<ByteType>();
		edgelPatchMasks = bfactory.create(img, new ByteType());

	}

	public ArrayList<Edgel> getEdgels(){
		return edgels;
	}
	
	public void genImgChunks() {

		if (imgChunks == null) {
			ImgFactory<B> factory = mask.factory();
			imgChunks = factory.create(mask, mask.firstElement());
		}
	}

	public void computeEdgels() {

		edgels = SubpixelEdgelDetection.getEdgels(mask, mask.factory(), 10.0f);
		logger.info("num edgels " + edgels.size());

	}
	
	public void setEdgels(ArrayList<Edgel> edgels){
		this.edgels = edgels;
	}

	public void crackBoundaries() {
		// compute crack orientation:
		ImgMoment imgmom = new ImgMoment();
		double[] or = imgmom.orientation(mask);
		imgmom.orientationEvalsEvecs(or);

		double[] cent = imgmom.getM1();
		ArrayUtil.divide(cent, imgmom.getM0());

		//System.out.println(" pos: " + ArrayUtil.printArray( ArrayUtil.reshape2D( or , 3, 3, false )));
		logger.debug(" cent:\n " + ArrayUtil.printArray(cent));
		logger.debug(" orEvecs:\n "
				+ ArrayUtil.printArray(ArrayUtil.reshape2D(imgmom.getEvecs(),
						3, 3, false)));

	}

	/**
	 * Returns the index k into the edgel list closest to the input
	 * point.
	 */
	public int edgelIdxNearest(float[] pos) {
		int i = -1;

		float minDist = Float.MAX_VALUE;
		int n = 0;
		for (Edgel e : edgels) {
			float f = ArrayUtil.sumSquares(ArrayUtil.subtract(pos,
					e.getPosition()));
			if (f < minDist) {
				i = n;
				minDist = f;
			}
			n++;
		}

		return i;
	}
	
	/**
	 * Returns the integer coordinate of the midpoint of a patch
	 * with the specified size.  Is only accurate for patches
	 * with an odd size in every dimension.
	 * 
	 * @param patchSize size of the patch
	 * @return the midpoint coordinate
	 */
	public static int[] patchSizeToMidpt(int[] patchSize){ // determine translation

		int[] midPt = ArrayUtil.clone(patchSize);
		ArrayUtil.addInPlace(midPt, -1);
		ArrayUtil.divide(midPt, 2);

		return midPt;
	}
	
	public Img<T> edgelToImage(Edgel edgel, Img<T> src, int[] patchSize) {

		logger.info(" edgel pos : " + ArrayUtil.printArray(edgel.getPosition()));
		logger.info(" edgel grad: " + ArrayUtil.printArray(edgel.getGradient()));
		logger.info(" edgel mag : " + edgel.getMagnitude());

		AffineTransform3D xfm = pickTransformation(edgel);
		Img<T> res = normalPatch(ArrayUtil.toDouble(edgel.getPosition()),
				patchSize, xfm, src);

		CrackCorrection.setMask(ArrayUtil.toDouble(edgel.getPosition()),
				patchSize, xfm, edgelPatchMasks, new ByteType(b));
		b++;

		ImgUtil.printNumNonZero(this.edgelPatchMasks);
		//		System.out.println("nnz edgelPatch: " + );

		return res;
	}
	
	public <M extends RealType<M>> void maskToDepthSample( RandomAccessible<M> maskView, RandomAccessible<T> imgView, int[] patchSize, 
			double thresh, boolean greaterThanThresh  )
	{
		
		RandomAccess<M> funcra = maskView.randomAccess();
		RandomAccess<T> destra = imgView.randomAccess();
		
		int ndims 	  = img.numDimensions();
		int ndims_out = patchSize.length;
		int zdimIdx = ndims - 1;

		int[] patchSizeAug = new int[ndims];
		for (int i=0; i<ndims_out; i++){
			patchSizeAug[i]=patchSize[i];
		}
		if(ndims!=ndims_out){
			for (int i=ndims_out; i<ndims; i++){
				patchSizeAug[i] = 1;
			}
		}

		int[] midPt = patchSizeToMidpt(patchSizeAug);
		funcra.setPosition(midPt);

		double[] testpos = new double[3];
		funcra.localize(testpos);
		logger.debug(" edgel patch position: " + ArrayUtil.printArray(testpos));


		double ctrVal = (funcra.get()).getRealDouble();

		funcra.fwd(zdimIdx);
		double fwdVal = (funcra.get()).getRealDouble();

		funcra.localize(testpos);
		logger.debug(" edgel patch position: " + ArrayUtil.printArray(testpos));


		funcra.bck(zdimIdx);
		funcra.bck(zdimIdx);
		double bckVal = (funcra.get()).getRealDouble();

		funcra.localize(testpos);
		logger.trace(" edgel patch position: " + ArrayUtil.printArray(testpos));

		logger.trace(" dist val at bck: " + bckVal );
		logger.trace(" dist val at ctr: " + ctrVal );
		logger.trace(" dist val at fwd: " + fwdVal );

		boolean goFwd = false;
		if( fwdVal < bckVal ){
			logger.info( " GOOD " );
		}else if( fwdVal > ctrVal ){
			logger.info( " BAD " );
			goFwd = true;
		}else{
			//they're equal
			logger.warn( "CAN NOT DETERMINE DIRECTION" );
		}

		int[] debugPatchPos = new int[imgPatch.numDimensions()]; 
		int[] debugDra = new int[funcra.numDimensions()];

		double d = -1;
		double dist = 1;
		Cursor<T> curs          = imgPatch.cursor();
		Cursor<FloatType> dcurs = depthPatch.cursor();
		while(curs.hasNext())
		{
			curs.fwd();
			dcurs.fwd();

			curs.localize(debugPatchPos);

			logger.trace(" patch pos: " + ArrayUtil.printArray(debugPatchPos));

			funcra.setPosition(curs);
			funcra.setPosition(0, zdimIdx);

			funcra.localize(debugDra);
			logger.trace(" distra pos : " + ArrayUtil.printArray(debugDra));


			dist = 1;
			for(int n=0; n<maxNumNormalSteps; n++){

				if(goFwd){
					funcra.fwd(zdimIdx);
				}else{
					funcra.bck(zdimIdx);
				}

				funcra.localize(debugDra);
				logger.trace(" distra pos : " + ArrayUtil.printArray(debugDra));

				d = (funcra.get()).getRealDouble(); 


				logger.trace(" func val here " + d );

				if(greaterThanThresh){

					if( d < thresh){ // we're outside the crack
						dist -= d;
						break;
					}else{
						dist ++;
					}
					
				}else{
					
					if( d > thresh){ // we're outside the crack
						dist -= d;
						break;
					}else{
						dist ++;
					}
				}
			}

			// set depth
			dcurs.get().setReal(dist);

			// set image
			destra.setPosition(funcra);
			curs.get().setReal( (destra.get()).getRealDouble() );
			

//			break; // for debug only 
			
		}

		logger.trace(" dist comp w normality: " + dist );
	}

	public void computeCrackDepthNormal(Edgel edgel, int[] patchSize)
	{
	
		 //edgelToImage(edgel, img, patchSize);
		 
		 		
		 int ndims 	  = img.numDimensions();
		int ndims_out = patchSize.length;
		
		int[] patchSizeAug = new int[ndims];
		for (int i=0; i<ndims_out; i++){
			patchSizeAug[i]=patchSize[i];
		}
		if(ndims!=ndims_out){
			for (int i=ndims_out; i<ndims; i++){
				patchSizeAug[i] = 1;
			}
		}
		
		if( imgPatch == null){
			imgPatch = img.factory().create(patchSize, img.firstElement());
		}
		if( depthPatch == null){
			ArrayImgFactory<FloatType> ffactory = new ArrayImgFactory<FloatType>();
			depthPatch = ffactory.create(patchSize, new FloatType());
		}
		
		
//		int[] patchSizeWChan = new int[patchSize.length + 1];
//		for(int i=0; i<patchSize.length; i++)
//		{
//			patchSizeWChan[i] = patchSize[i];
//		}
//		patchSizeWChan[patchSize.length] = nChannels;
		
		
		RealTransformRandomAccessible<T,?> imgEdgelView = edgelToView( edgel, img, patchSizeAug );
		RealTransformRandomAccessible<B,?> mskEdgelView = edgelToView( edgel, mask, patchSizeAug );
//		RandomAccess<T> ra = rtra.randomAccess();

//		Img<FloatType> maskDist = ImgUtil.signedDistance(mask);
//		
//		RealTransformRandomAccessible<FloatType,?> distrtra = edgelToView( edgel, maskDist, patchSizeAug );
//		RandomAccess<FloatType> distra = distrtra.randomAccess();
//		distToMaskSDF( distra, ra, patchSize, 0, true );
		
		
		maskToDepthSample( mskEdgelView, imgEdgelView, patchSize, 128, false );
		
		
	}
	
	public Img<FloatType> compMaskDistMgdm()
	{
		int[][][] maskInt = ImgUtil.toIntArray3d(mask);
		boolean[][][] m = new boolean[maskInt.length][maskInt[0].length][maskInt[0][0].length];
		ArrayUtil.fill(m, true);

		MgdmDecomposition mgdm = new MgdmDecomposition(maskInt, 2, m);
		float[][][] d1 = mgdm.exportDistanceForObject3d(0);

		Img<FloatType> distXfmImg = null;

		//			logger.info( "writing patches" );
		ArrayImgFactory<FloatType> factory = new ArrayImgFactory<FloatType>();

		distXfmImg = factory.create(img, new FloatType(img
				.firstElement().getRealFloat()));
		ImgUtil.copyToImg(distXfmImg, d1);
		
		ImagePlus ip1;
//		try {
//			
//			ip1 = ImgUtil.toImagePlus(distXfmImg);
//			
//			IJ.save(ip1,
//					"/groups/jain/home/bogovicj/projects/crackSegmentation/toy/distXfm.tif");
//
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
		
		return distXfmImg;
	}

	public Img<T> computeLocalCrackDepthDistFunc(Edgel edgel, int[] patchSize) {
		logger.debug("computeLocalCrackDepthDistFunc");

//		ImgUtil.printNumNonZero(mask);
//		ArrayList<Integer> uniqueList = ImgUtil.uniqueInt(mask);
//		System.out.println("unique in mask:\n" + uniqueList);

		int[][][] maskInt = ImgUtil.toIntArray3d(mask);
		boolean[][][] m = new boolean[maskInt.length][maskInt[0].length][maskInt[0][0].length];
		ArrayUtil.fill(m, true);

		MgdmDecomposition mgdm = new MgdmDecomposition(maskInt, 2, m);
		float[][][] d1 = mgdm.exportDistanceForObject3d(0);

		try {

			//			logger.info( "writing patches" );
			ArrayImgFactory<FloatType> factory = new ArrayImgFactory<FloatType>();

			Img<FloatType> distXfmImg = factory.create(img, new FloatType(img
					.firstElement().getRealFloat()));

			ImgUtil.copyToImg(distXfmImg, d1);
			ImagePlus ip1 = ImgUtil.toImagePlus(distXfmImg);
			IJ.save(ip1,
					"/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/distXfm.tif");

		} catch (Exception e) {
			e.printStackTrace();
		}

//		int c =  ImgUtil.numLessThanZero(d1);
//		System.out.println(" count dist<=0 :" + c);

		Img<T> d1img = img.factory().create(img, img.firstElement());
		ImgUtil.copyToImg(d1img, d1);

		Img<T> depthPatch = edgelToImage(edgel, d1img, patchSize);

//		int cp = ImgUtil.numLessThanZero(depthPatch);
//		System.out.println(" count patch dist<=0 :" + cp);

		return depthPatch;
	}



	public static void testCrackCorrPipeline() {

		logger.debug("testCrackCorrPipeline");

		String imgfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/crackVolDown_cp.tif";
		
//		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp.tif";
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_smooth.tif";
//		String mgdmfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/LabelsMgdm_ds_interp_cp_v2.tif";

		//		String maskWSfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/Labels_ds_interp_cp_manPaint2.tif";
		String maskWSfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/Labels_ds_interp_cp_manPaint3.tif";
		
		String depthPatchFn = "/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/depthPatch_nrm";
		String imgPatchFn = "/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/imgPatch_nrm";
		

		String patchOut = "/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/patch";
		String edgelPatchOut = "/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/edgel_patch";
		String depthPatchOut = "/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/edgel_depthpatch";
		
		String edgelMaskPrefix = "/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/edgelMask";

		//		ImagePlus imgip  = IJ.openImage(imgfn);
		//		ImagePlus maskip = IJ.openImage(maskWSfn);

		//		Img<FloatType> img = ImagePlusAdapter.convertFloat(imgip);
		//		ByteImagePlus<UnsignedByteType> mask = ImagePlusAdapter.wrapByte(maskip);

		ArrayImgFactory<UnsignedByteType> ubfactory = new ArrayImgFactory<UnsignedByteType>();
		ArrayImgFactory<FloatType> ffactory = new ArrayImgFactory<FloatType>();

		Img<FloatType> img = null;
		Img<UnsignedByteType> mask = null;
		try {
			img = new ImgOpener().openImg(imgfn, ffactory, new FloatType());
			mask = new ImgOpener().openImg(maskfn, ubfactory,
					new UnsignedByteType());
		} catch (ImgIOException e) {
			e.printStackTrace();
		}

		CrackCorrection<FloatType, UnsignedByteType> cc = new CrackCorrection<FloatType, UnsignedByteType>(
				img, mask);

		cc.computeEdgels();

		int[] N = new int[] { 19, 19, 7 };

//		float[] tgtPos  = new float[]{ 66f, 290f, 13f };
//		float[] tgtPos = new float[]{ 70f, 312f, 13f };
//		float[] tgtPos = new float[] { 139f, 227f, 13f };
		float[] tgtPos = new float[] { 142, 221f, 16f };

		int i = cc.edgelIdxNearest(tgtPos);
		String fnSuffix = "_"+i+".tif" ;
		
		logger.info("\nedgel idx: " + i);

//		Img<FloatType> img1 = cc.edgelToImage(cc.edgels.get(i), img, N);
//		Img<FloatType> depthPatch = cc.computeLocalCrackDepthDistFunc(
//				cc.edgels.get(i), N);
//		ip1 = ImgUtil.toImagePlus(img1);
//		IJ.save(ip1, patchOut + "_" + i + ".tif");
//		ImagePlus ipdp = ImgUtil.toImagePlus( depthPatch );
//		IJ.save(ipdp, depthPatchOut + "_" + i + ".tif");
		
		cc.computeCrackDepthNormal( cc.edgels.get(i), N );
		
		ImgUtil.write( cc.depthPatch, 		depthPatchFn + fnSuffix );
		ImgUtil.write( cc.imgPatch, 		imgPatchFn + fnSuffix );
		ImgUtil.write( cc.edgelPatchMasks, 	edgelMaskPrefix + fnSuffix );
		
		

//		try 
//		{
//			ImagePlus ipep = ImgUtil.toImagePlus( cc.edgelPatchMasks );
//			IJ.save(ipep, edgelPatchOut + "_" + i + ".tif");
//
//		} 
//		catch (Exception e) 
//		{
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
	}

	public static void testDistanceBasedCrackSideComp() {

		String imgfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/crackVolDown_cp.tif";
		//		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp.tif";
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/Labels_ds_interp_cp_manPaint3.tif";

		String mgdmfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/LabelsMgdm_ds_interp_cp_v2.tif";
		String mgdmdstfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/DistMgdm_ds_interp_cp_v2.tif";

		String maskWSfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_manPaint2.tif";

		ImagePlus imgip = IJ.openImage(imgfn);
		ImagePlus maskip = IJ.openImage(maskfn);

		//		System.out.println("imgip " + imgip);
		//		System.out.println("maskip" + maskip);

		Img<FloatType> img = ImagePlusAdapter.convertFloat(imgip);
		ByteImagePlus<UnsignedByteType> mask = ImagePlusAdapter
				.wrapByte(maskip);

		//		System.out.println("img " + img);
		//		System.out.println("mask " + mask);

		CrackCorrection<FloatType, UnsignedByteType> cc = new CrackCorrection<FloatType, UnsignedByteType>(
				img, mask);

		//cc.computeEdgels();
		//cc.crackBoundaries();

		boolean[][][] mskArray = ImgUtil.toBooleanArray3dNeg(mask);
		//float[][][] imgArray = toFloatArray3d(img);

		System.out.println("nnz in mask: \n"
				+ ArrayUtil.nnz(ArrayUtil.reshape1D(mskArray, true)));

		int[][][] sideArray = new int[mskArray.length][mskArray[0].length][mskArray[0][0].length];

		sideArray[0][0][13] = 1;
		sideArray[287][437][13] = 2;

		/***********/
		//		MgdmDecomposition mgdm = new MgdmDecomposition(sideArray,3, mskArray);
		//
		//		int[][][][] mgdmLabels = mgdm.exportAllLabels3d();
		//
		//		Img<UnsignedByteType> mgdmLImg = mask.factory().create(
		//				new long[]{ mask.dimension(0), mask.dimension(1),
		//					mask.dimension(2), 3 },
		//					new UnsignedByteType());

		//	    Img<FloatType> d1img = img.factory().create(img, img.firstElement());
		//	    Img<FloatType> d2img = img.factory().create(img, img.firstElement());
		//
		//	    float[][][] d1 = mgdm.exportDistanceForObject3d(1);
		//	    float[][][] d2 = mgdm.exportDistanceForObject3d(2);

		try {

			//			ImgUtil.copyToImg4d(mgdmLImg, mgdmLabels);
			//			ImagePlus ipMgdm = ImgUtil.toImagePlus(mgdmLImg);
			//			IJ.save(ipMgdm, mgdmfn);

			//	       copyToImg4d(d1img, d1);
			//	       ImagePlus ipMgdm = toImagePlus(d1img);
			//	       IJ.save(ipMgdm, d1fn);

			//				       copyToImg4d(d2img, d2);
			//				       ipMgdm = toImagePlus(d2img);
			//				       IJ.save(ipMgdm, d2fn);

			//			ImgUtil.copyToImg( d1img, d1);
			//			ImagePlus ip = ImgUtil.toImagePlus(d1img);
			//			IJ.save( ip, mgdmdstfn );

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static void testEdgelResamp2() {

		System.out.println("test Edgel Resamp");

		//Img<FloatType> testImg = ImgUtil.createGradientImgX( 20, 20, 20, new FloatType() );
		//Img<FloatType> testImg = ImgUtil.createGradientImgY( 20, 20, 20, new FloatType() );
		Img<FloatType> testImg = ImgUtil.createGradientImgZ(20, 20, 20,
				new FloatType());

		CrackCorrection<FloatType, UnsignedByteType> cc = new CrackCorrection<FloatType, UnsignedByteType>(
				testImg, null);

		//float[] norm = new float[]{ 1.0f, 1.0f, 0.0f };
		float[] norm = new float[] { 0.0f, 0.0f, 1.0f };
		norm = ArrayUtil.normalizeLength(norm);

		System.out.println(" normal vector: (" + ArrayUtil.printArray(norm)
				+ ") ");

		Edgel edgel = new Edgel(new float[] { 10.5f, 10.5f, 10.5f }, norm, 1);

		int[] patchSize = new int[] { 5, 5, 3 };
		Img<FloatType> resImg = cc.edgelToImage(edgel, testImg, patchSize);

		//ImageJFunctions.show(resImg);

	}

	public static void testEdgelResamp() {

		System.out.println("test Edgel Resamp");

		//Img<FloatType> testImg = ImgUtil.createGradientImgX( 20, 20, 20, new FloatType() );
		//Img<FloatType> testImg = ImgUtil.createGradientImgY( 20, 20, 20, new FloatType() );
		Img<FloatType> testImg = ImgUtil.createGradientImgZ(20, 20, 20,
				new FloatType());

		System.out.println(" testImg: " + testImg);

		//float[] norm = new float[]{ 1.0f, 1.0f, 0.0f };
		float[] norm = new float[] { 0.0f, 0.0f, 1.0f };
		norm = ArrayUtil.normalizeLength(norm);

		System.out.println(" normal vector: (" + ArrayUtil.printArray(norm)
				+ ") ");

		Edgel edgel = new Edgel(new float[] { 10.5f, 10.5f, 10.5f }, norm, 1);
		System.out.println("edgel : " + edgel);

		NLinearInterpolatorFactory<FloatType> interpFactory = new NLinearInterpolatorFactory<FloatType>();

		RealRandomAccessible<FloatType> interpolant = Views.interpolate(
				Views.extendMirrorSingle(testImg), interpFactory);

		RealRandomAccess<FloatType> rra = interpolant.realRandomAccess();
		double[] pos = new double[] { 10.5, 10.5, 10.5 };
		rra.setPosition(pos);
		double val = rra.get().getRealDouble();

		System.out.println(" value at: (" + ArrayUtil.printArray(pos) + "): "
				+ val);

		AffineTransform3D xfm = pickTransformation(edgel);

		int[] N = new int[] { 5, 5, 3 };
		Img<FloatType> resImg = normalPatch(pos, N, xfm, testImg);

		//ImageJFunctions.show(resImg);

	}

	public static AffineTransform3D edgelToXfm(Edgel edgel, int[] midPt)
	{
		AffineTransform3D rotXfm = CrackCorrection.pickTransformation(edgel);
		
		/******************************/
		AffineTransform3D xfm = rotXfm.copy();

		int ndims_in = midPt.length;
		double[] target = new double[ndims_in];


		xfm.apply( ArrayUtil.toDouble(midPt), target );
		double[] diff = ArrayUtil.subtract(
				ArrayUtil.toDouble(edgel.getPosition()), target);

		for (int i = 0; i < ndims_in; i++) {
			xfm.set(diff[i], i, ndims_in);
		}
		/******************************/
		
		return xfm;	
	}

	
	/**
	 * Returns a transformed view of the input source image relative to an edgel.
	 * The {@link Edgel} position will map to the midpoint the output patch.
	 * The +z axis of the output view corresponds to the gradient direction of 
	 * the edgel.
	 * 
	 * @param edgel the edgel
	 * @param src the source image
	 * @param patchSize the patch size
	 * @return the transformed view into the source image
	 */
	public static <T extends RealType<T>> RealTransformRandomAccessible<T,InverseRealTransform> edgelToView(Edgel edgel, Img<T> src, int[] patchSize) 
	{
		logger.info(" edgel pos : " + ArrayUtil.printArray(edgel.getPosition()));
		logger.info(" edgel grad: " + ArrayUtil.printArray(edgel.getGradient()));
		logger.info(" edgel mag : " + edgel.getMagnitude());
		
		int[] midPt = patchSizeToMidpt( patchSize ); 	
		AffineTransform3D xfm = edgelToXfm(edgel, midPt);


		NLinearInterpolatorFactory<T> interpFactory = new NLinearInterpolatorFactory<T>();

		RealRandomAccessible<T> interpolant = Views.interpolate(
				Views.extendMirrorSingle(src), interpFactory);

		RealTransformRandomAccessible<T, InverseRealTransform> rv = 
				RealViews.transform( interpolant, xfm.inverse() );
		
		return rv;
	}
	
	// Views.offset( img, 1, 0 0 )
	
	/**
	 *  Samples a patch of size N^d where d is the dimensionality of the source image src. 
	 *  The input transformation maps the x- and y-axes to the axes of the sampling plane.
	 *  It is typically determined by the {@link pickTransformation} method.
	 *
	 * @param basepos point
	 * @param N size of a dimension of the output patch
	 * @param xfm transformation
	 * @param src source image
	 * @return
	 */
	public static <T extends RealType<T>> Img<T> normalPatch(double[] basepos,
			int[] patchSize, AffineTransform3D xfmIn, Img<T> src) {

		AffineTransform3D xfm = xfmIn.copy();

		NLinearInterpolatorFactory<T> interpFactory = new NLinearInterpolatorFactory<T>();

		RealRandomAccessible<T> interpolant1 = Views.interpolate(
				Views.extendMirrorSingle(src), interpFactory);
		RealRandomAccess<T> srcRealRa = interpolant1.realRandomAccess();

		int ndims_in = src.numDimensions();

		Img<T> out = src.factory().create(patchSize,
				Util.getTypeFromRandomAccess(src));

		logger.debug(" basepos  :" + ArrayUtil.printArray(basepos));
		logger.debug(" out size :" + ArrayUtil.printArray(patchSize));

		int[] pos = new int[ndims_in];
		double[] dpos = new double[ndims_in];

		//		double[][][] a = new double[patchSize[0]][patchSize[1]][patchSize[2]];

		// determine translation
		double[] midPt = ArrayUtil.toDouble(patchSize);
		ArrayUtil.addInPlace(midPt, -1);
		ArrayUtil.divide(midPt, 2);
		logger.debug(" midPt :" + ArrayUtil.printArray(midPt));

		double[] target = new double[ndims_in];

		xfm.apply(midPt, target);
		double[] diff = ArrayUtil.subtract(basepos, target);

		logger.debug(" diff  :" + ArrayUtil.printArray(diff));

		for (int i = 0; i < ndims_in; i++) {
			xfm.set(diff[i], i, ndims_in);
		}

		//		double[][] c = new double[ndims_in][ndims_in+1];
		//		xfm.toMatrix(c);
		//		logger.debug("xfm: " + ArrayUtil.printArray(c));

		Cursor<T> cursor = out.cursor();
		while (cursor.hasNext()) {

			cursor.fwd();
			cursor.localize(pos);

			double[] pt = ArrayUtil.toDouble(pos);

			xfm.apply(pt, dpos);

			logger.trace(" pt     :" + ArrayUtil.printArray(pt));
			logger.trace(" dpos   :" + ArrayUtil.printArray(dpos));

			srcRealRa.setPosition(dpos);
			double val = srcRealRa.get().getRealDouble();
			cursor.get().setReal(val);

			logger.trace(" val    :" + val);


		}

		//		System.out.println("res:\n" + ArrayUtil.printArray(a));

		return out;
	}

//	/**
//	 *  Samples a patch of size N^d where d is the dimensionality of the
//	 *  source image src. The input transformation maps the x- and y-axes
//	 *  to the axes of the sampling plane and is typically determined by
//	 *  the {@link pickTransformation} method.
//	 *
//	 *
//	 * @param basepos point
//	 * @param N size of a dimension of the output patch
//	 * @param xfm transformation
//	 * @param src source image
//	 * @return
//	 */
//	public static <T extends RealType<T>> Img<T> normalPatchCoords(
//			double[] basepos, int[] patchSize, AffineTransform3D xfmIn,
//			Img<T> src) {
//
//		int ndims_in = src.numDimensions();
//
//		AffineTransform3D xfm = xfmIn.copy();
//		int[] imgOutSz = new int[patchSize.length + 1];
//
//		for (int i = 0; i < patchSize.length; i++) {
//			imgOutSz[i] = patchSize[i];
//		}
//		imgOutSz[patchSize.length] = ndims_in;
//
//		Img<T> out = src.factory().create(patchSize,
//				Util.getTypeFromRandomAccess(src));
//		CompositeIntervalView<T, ?> outCol = Views.collapse(out);
//
//		int[] 	 pos  = new int   [ndims_in];
//		double[] dpos	= new double[ndims_in];
//		
//		// determine translation
//		double[] midPt = ArrayUtil.toDouble(patchSize);
//		ArrayUtil.addInPlace(midPt, -1);
//		ArrayUtil.divide(midPt, 2);
//		logger.debug(" midPt :" + ArrayUtil.printArray( midPt ) );
//		
//		double[] target = new double[ndims_in];
//		
//		xfm.apply(midPt, target);
//		double[] diff = ArrayUtil.subtract( basepos, target );
//		
//		logger.debug(" diff  :" + ArrayUtil.printArray( diff ) );
//		
//		for( int i=0; i<ndims_in; i++ ){
//			xfm.set( diff[i], i, ndims_in);
//		}
//		
//		Cursor<T> cursor = out.cursor();
//		Cursor<T> outCol = outCol.cursor();
//		while(cursor.hasNext()){
//			
//			cursor.fwd();
//			cursor.localize(pos);
//			//cursor.setP
//			
//			double[] pt = ArrayUtil.toDouble(pos);
//			
//			xfm.apply( pt, dpos);
//			
//			logger.debug(" pt     :" + ArrayUtil.printArray( pt  ) );
//			logger.debug(" dpos   :" + ArrayUtil.printArray( dpos) );
//			
//			cursor.get().setReal( val );
//			
//			logger.debug(" val    :" + val );
//			
//		}
//		
//		return out;
//	}
	
	
	/**
	 *  Samples a patch of size N^d where d is the dimensionality of the 
	 *  source image src. The input transformation maps the x- and y-axes 
	 *  to the axes of the sampling plane and is typically determined by 
	 *  the {@link pickTransformation} method.
	 * 
	 *
	 * @param basepos point
	 * @param N size of a dimension of the output patch
	 * @param xfm transformation
	 * @param src source image
	 * @return 
	 */
	public static <L extends NumericType<L>> void setMaskXfm(double[] basepos, int[] patchSize,  
			AffineTransform3D xfmIn,  RandomAccessible<L> src, L val ){
		
		AffineTransform3D xfm = xfmIn.copy();
		
		RandomAccess<L> ra = src.randomAccess();
//		ra.setPosition(position);
		
		int ndims_in  = src.numDimensions();
		
		logger.debug(" basepos  :" + ArrayUtil.printArray( basepos ) );
		logger.debug(" out size :" + ArrayUtil.printArray( patchSize ) );

		double[]  pos	= new double[ndims_in];
		double[] dpos	= new double[ndims_in];
		
		// determine translation
		double[] midPt = ArrayUtil.toDouble(patchSize);
		ArrayUtil.addInPlace(midPt, -1);
		ArrayUtil.divide(midPt, 2);
		logger.debug(" midPt :" + ArrayUtil.printArray( midPt ) );
		
		double[] target = new double[ndims_in];
		
		xfm.apply(midPt, target);
		double[] diff = ArrayUtil.subtract( basepos, target );
		
		logger.debug(" diff  :" + ArrayUtil.printArray( diff ) );
		
		for( int i=0; i<ndims_in; i++ ){
			xfm.set( diff[i], i, ndims_in);
		}

		IntervalIterator iter = new IntervalIterator(patchSize);
		while(iter.hasNext()){
			
			iter.fwd();
			iter.localize(pos);
			
			xfm.apply( pos, dpos);
			int[] pt = ArrayUtil.toIntRound(dpos);
			logger.trace(" pt :" + ArrayUtil.printArray( pt ) );
			
			ra.setPosition(pt);
			ra.get().set(val);
		}
	}
	
	
	/**
	 *  Samples a patch of size N^d where d is the dimensionality of the 
	 *  source image src. The input transformation maps the x- and y-axes 
	 *  to the axes of the sampling plane and is typically determined by 
	 *  the {@link pickTransformation} method.
	 * 
	 *
	 * @param basepos point
	 * @param N size of a dimension of the output patch
	 * @param xfm transformation
	 * @param src source image
	 * @return 
	 */
	public static <L extends NumericType<L>> void setMask(double[] basepos, int[] patchSize,  
			AffineTransform3D xfmIn,  RandomAccessible<L> src, L val ){
		
		AffineTransform3D xfm = xfmIn.copy();
		
		RandomAccess<L> ra = src.randomAccess();
//		ra.setPosition(position);
		
		int ndims_in  = src.numDimensions();
		
		logger.debug(" basepos  :" + ArrayUtil.printArray( basepos ) );
		logger.debug(" out size :" + ArrayUtil.printArray( patchSize ) );

		double[]  pos	= new double[ndims_in];
		double[] dpos	= new double[ndims_in];
		

		IntervalIterator iter = new IntervalIterator(patchSize);
		while(iter.hasNext()){
			
			iter.fwd();
			iter.localize(pos);
			
			xfm.apply( pos, dpos);
			int[] pt = ArrayUtil.toIntRound(dpos);
			logger.trace(" pt :" + ArrayUtil.printArray( pt ) );
			
			ra.setPosition(pt);
			ra.get().set(val);
		}
	}
	
	public static AffineTransform3D pickTransformation(Edgel edgel){
	
		float[] nrm = edgel.getGradient();
		
		DenseMatrix64F normMtx = new DenseMatrix64F( new double[][]{ ArrayUtil.toDouble(nrm) } );
		DenseMatrix64F remMtx = remainderSubspace(normMtx);
		logger.debug("rem Mtx: " + remMtx );
		
		double[][] Ra = new double[3][4];
		for(int i=0; i<3; i++){
			for(int j=0; j<2; j++){
				Ra[i][j] = remMtx.get( i, j );
			}
			Ra[i][2] = normMtx.get(i);
		}
		
		logger.info(" Ra: \n" + ArrayUtil.printArray(Ra) );
		
		AffineTransform3D R = new AffineTransform3D();
		R.set(Ra);
		
//		AffineTransform3D R = RodriguesRotation.rotation(
//		new double[]{0,0,1}, ArrayUtil.toDouble(nrm));
		
		return R;
	}

   /**
	 * Replace this with a cross product if too slow..
	 * @param in an M x N matrix whose rows span a subspace of R^N
	 * @return an orthogonal basis for the subspace of R^N not spanned by in
	 */
	public static DenseMatrix64F remainderSubspace(DenseMatrix64F in){
		SimpleSVD<SimpleMatrix> svd = new SimpleSVD<SimpleMatrix>(in,false);
		SimpleMatrix nullspace = svd.nullSpace();
		return nullspace.getMatrix();
	}
	

	
	public static void testCombReplace(){
		
//		HashSet<UnsignedByteType> set = new HashSet<UnsignedByteType>();
//		set.add(new UnsignedByteType((byte)0));
//		set.add(new UnsignedByteType((byte)1));
//		
//		System.out.println("unique values:\n" + set);
		
//		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp.tif";
//		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/Labels_ds_interp_cp_manPaint2.tif";
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/intermRes/Labels_ds_interp_cp_manPaint3.tif";

		IntType type = new IntType();
		ArrayImgFactory< IntType > factory = new ArrayImgFactory< IntType >();
		
		Img<IntType> mask = null;
		try 
		{
			mask = new ImgOpener().openImg( maskfn , factory, type );
		}
		catch (ImgIOException e)
		{
			e.printStackTrace();
		}
		
		ImgUtil.printNumNonZero(mask);
		
		ArrayList<Integer> uniqueSet = ImgUtil.uniqueInt(mask);
		System.out.println("unique values: " + uniqueSet.size());
		System.out.println("unique values:\n" + uniqueSet);
		
		
	}
	
	   
	public static void main(String[] args){

//		testDistanceBasedCrackSideComp();	
		
//		testEdgelResamp();
		
//		testEdgelResamp2();
		
		testCrackCorrPipeline();
		
//		testCombReplace();
		
		System.out.println("crack correction finished");
		System.exit(0);
	}


}
