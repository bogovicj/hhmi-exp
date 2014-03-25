package net.imglib2.algorithms.crack;

import java.util.ArrayList;
import java.util.Arrays;

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
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.img.imageplus.*;
import net.imglib2.algorithm.edge.*;
import net.imglib2.algorithms.geometry.RodriguesRotation;
import net.imglib2.algorithms.moments.ImgMoment;
import net.imglib2.img.*;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.InvertibleRealTransform;
import net.imglib2.type.numeric.*;
import net.imglib2.type.numeric.integer.*;
import net.imglib2.type.numeric.real.*;
import net.imglib2.view.Views;
import net.imglib2.util.ImgUtil;
import net.imglib2.util.Util;


/**
 * 
 * Class performing correction of high-pressure freezing cracks 
 * in electron microscopy images.
 * 
 * @author John Bogovic
 *
 * @param <T> image type
 * @param <B> crack label type 
 */
public class CrackCorrection<T extends RealType<T>, B extends AbstractIntegerType<B>> {

	Img<T> img;
	Img<B> mask;
	
	RealRandomAccess<T> interpRa;

	//IntegralImgDouble<B> maskIntImg;
	ArrayList<Edgel> edgels;
	
	Img<B> imgChunks;
	
	protected static Logger logger = LogManager.getLogger(CrackCorrection.class.getName());

	public CrackCorrection(){ }

	public CrackCorrection( Img<T> img, Img<B> mask ){
		this.img = img;
		this.mask = mask;

		NLinearInterpolatorFactory<T> interpFactory = 
				new NLinearInterpolatorFactory<T>();

		RealRandomAccessible<T> interpolant = 
				Views.interpolate( 
						Views.extendMirrorSingle( img ), interpFactory
				);
		
		interpRa = interpolant.realRandomAccess();
	}

	public void genImgChunks(){

		if(imgChunks == null){
			ImgFactory<B> factory = mask.factory();	
			imgChunks = factory.create(mask, mask.firstElement());
		}
	}

	public void computeEdgels(){
		
		edgels = SubpixelEdgelDetection.getEdgels(mask, mask.factory(), 2f);
		logger.info("num edgels " + edgels.size());
		
	}

	public void crackBoundaries(){
		// compute crack orientation:
		ImgMoment imgmom = new ImgMoment();
		double[] or = imgmom.orientation( mask );	
		imgmom.orientationEvalsEvecs( or );

		double[] cent = imgmom.getM1();
		ArrayUtil.divide(cent, imgmom.getM0());

		//System.out.println(" pos: " + ArrayUtil.printArray( ArrayUtil.reshape2D( or , 3, 3, false )));
		logger.debug(" cent:\n " + ArrayUtil.printArray( cent ));
		logger.debug(" orEvecs:\n " + ArrayUtil.printArray( ArrayUtil.reshape2D( imgmom.getEvecs(), 3, 3, false )));	


	}

	/**
	 * Returns the index k into the edgel list closest to the input
	 * point.
	 */
	public int edgelIdxNearest(float[] pos){
		int i = -1;

		float minDist = Float.MAX_VALUE;
		int n = 0;
		for( Edgel e : edgels ){
			float f = ArrayUtil.sumSquares( ArrayUtil.subtract(pos, e.getPosition()));
			if( f < minDist ){
				i = n;
				minDist = f;
			}
			n++;
		}

		return i;
	}
	
	public Img<T> edgelToImage( Edgel edgel, Img<T> src, int[] patchSize ){
		
		AffineTransform3D xfm = pickTransformation(edgel);
		Img<T> res = normalPatch(
						ArrayUtil.toDouble(edgel.getPosition()), 
						patchSize, xfm, src 
					 );
		
		return res;
	}

	public void computeLocalCrackDepth(){
		logger.debug("computeLocalCrackDepth");	

	}
	

	public static void testCrackCorrPipeline(){

		logger.debug("testCrackCorrPipeline");

		String imgfn  = "/groups/jain/home/bogovicj/projects/crackSegmentation/crackVolDown_cp.tif";
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp.tif";
		String mgdmfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/LabelsMgdm_ds_interp_cp_v2.tif";
		
		String maskWSfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_manPaint2.tif";

		String patchOut1 = "/groups/jain/home/bogovicj/projects/crackSegmentation/patch1.tif";
		String patchOut2 = "/groups/jain/home/bogovicj/projects/crackSegmentation/patch2.tif";
		
		ImagePlus imgip = IJ.openImage(imgfn);
		ImagePlus maskip = IJ.openImage(maskWSfn);

		Img<FloatType> img = ImagePlusAdapter.convertFloat(imgip);
		ByteImagePlus<UnsignedByteType> mask = ImagePlusAdapter.wrapByte(maskip); 

		CrackCorrection<FloatType, UnsignedByteType> cc = 
				new CrackCorrection<FloatType, UnsignedByteType>( 
						img, mask );

		cc.computeEdgels();

		float[] tgtPos  = new float[]{ 66f, 290f, 13f };
		float[] tgtPos2 = new float[]{ 70f, 312f, 13f };

		int i = cc.edgelIdxNearest( tgtPos );
		int j = cc.edgelIdxNearest( tgtPos2 );

		logger.debug( "\n i: " + i );
		logger.debug( "ePos: " + ArrayUtil.printArray( cc.edgels.get(i).getPosition() ) );
		logger.debug( "eGrd: " + ArrayUtil.printArray( cc.edgels.get(i).getGradient() ) + "\n" );

		logger.debug( "\nj: " + j );
		logger.debug( "ePos: " + ArrayUtil.printArray( cc.edgels.get(j).getPosition() ) );
		logger.debug( "eGrd: " + ArrayUtil.printArray( cc.edgels.get(j).getGradient() ) + "\n" );


		AffineTransform3D xfm1 = pickTransformation(cc.edgels.get(i));
		AffineTransform3D xfm2 = pickTransformation(cc.edgels.get(j));

		int[] N = new int[]{ 19, 19, 7 };
		double[] p1 = ArrayUtil.toDouble(cc.edgels.get(i).getPosition());
		double[] p2 = ArrayUtil.toDouble(cc.edgels.get(j).getPosition());

		Img<FloatType> img1 = normalPatch( p1, N, xfm1, img);
	 	Img<FloatType> img2 = normalPatch( p2, N, xfm2, img);

		try {

			logger.debug( "writing patches" );
			ImagePlus ip1 = ImgUtil.toImagePlus(img1);
			IJ.save(ip1, patchOut1);

			ImagePlus ip2 = ImgUtil.toImagePlus(img2);
			IJ.save(ip2, patchOut2);

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static void testDistanceBasedCrackSideComp(){

		String imgfn  = "/groups/jain/home/bogovicj/projects/crackSegmentation/crackVolDown_cp.tif";
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp.tif";
		String mgdmfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/LabelsMgdm_ds_interp_cp_v2.tif";
		
		String maskWSfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_manPaint2.tif";

		ImagePlus imgip = IJ.openImage(imgfn);
		ImagePlus maskip = IJ.openImage(maskfn);

//		System.out.println("imgip " + imgip);
//		System.out.println("maskip" + maskip);

		Img<FloatType> img = ImagePlusAdapter.convertFloat(imgip);
		ByteImagePlus<UnsignedByteType> mask = ImagePlusAdapter.wrapByte(maskip); 

//		System.out.println("img " + img);
//		System.out.println("mask " + mask);

		CrackCorrection<FloatType, UnsignedByteType> cc = 
				new CrackCorrection<FloatType, UnsignedByteType>( 
						img, mask );

		//cc.computeEdgels();
		//cc.crackBoundaries();

		boolean[][][]   mskArray = ImgUtil.toBooleanArray3dNeg(mask);
		//float[][][] imgArray = toFloatArray3d(img);


		System.out.println("nnz in mask: \n" +
				ArrayUtil.nnz(ArrayUtil.reshape1D(mskArray, true))
				);

		int[][][]  sideArray = new int[mskArray.length][mskArray[0].length][mskArray[0][0].length];

		sideArray[0][0][13]     = 1;
		sideArray[287][437][13] = 2;

		MgdmDecomposition mgdm = new MgdmDecomposition(sideArray,3, mskArray);

		int[][][][] mgdmLabels = mgdm.exportAllLabels3d();

		Img<UnsignedByteType> mgdmLImg = mask.factory().create(
				new long[]{ mask.dimension(0), mask.dimension(1),
					mask.dimension(2), 3 },
					new UnsignedByteType());

//	    Img<FloatType> d1img = img.factory().create(img, img.firstElement());
//	    Img<FloatType> d2img = img.factory().create(img, img.firstElement());
//	   
//	    float[][][] d1 = mgdm.exportDistanceForObject3d(1);
//	    float[][][] d2 = mgdm.exportDistanceForObject3d(2);
	      
		try {

			ImgUtil.copyToImg4d(mgdmLImg, mgdmLabels);
			ImagePlus ipMgdm = ImgUtil.toImagePlus(mgdmLImg);
			IJ.save(ipMgdm, mgdmfn);

			//	       copyToImg4d(d1img, d1);
			//	       ImagePlus ipMgdm = toImagePlus(d1img);
			//	       IJ.save(ipMgdm, d1fn);
			//	       
			//	       copyToImg4d(d2img, d2);
			//	       ipMgdm = toImagePlus(d2img);
			//	       IJ.save(ipMgdm, d2fn);

		} catch (Exception e) {
			e.printStackTrace();
		}

	}


//	public static void testEdgelResamp(){
//
//		System.out.println("test Edgel Resamp");
//
//		Img<FloatType> testImg = createGradientImgX( 20, 20, 20, new FloatType() );
//		System.out.println(" testImg: " + testImg );
//	
//		float[] norm = new float[]{ 1.0f, 1.0f, 0.0f };
//		norm = ArrayUtil.normalizeLength(norm);
//
//		System.out.println(" normal vector: (" + ArrayUtil.printArray(norm) + ") "  );
//
//		Edgel edgel = new Edgel( new float[]{ 10.0f, 10.0f, 10.0f}, 
//								 norm, 1);
//		System.out.println("edgel : " + edgel);
//
//		NLinearInterpolatorFactory<FloatType> interpFactory = 
//			new NLinearInterpolatorFactory<FloatType>();
//
//		RealRandomAccessible<FloatType> interpolant = 
//			Views.interpolate( 
//					Views.extendMirrorSingle( testImg ), interpFactory
//			);
//
//		RealRandomAccess<FloatType> rra = interpolant.realRandomAccess();		
//		double[] pos = new double[]{ 10.5, 10.5, 10.5 } ;
//		rra.setPosition( pos );
//		double val = rra.get().getRealDouble();		
//
//		System.out.println(" value at: (" + ArrayUtil.printArray(pos) + "): " + val );
//
//		InvertibleRealTransform xfm = pickTransformation(edgel);
//		
//		RealTransformRandomAccessible<FloatType, InverseRealTransform> res = 
//				RealViews.transform(interpolant, xfm);
//		
//		RealTransformRealRandomAccess rtra = res.realRandomAccess();
//		rtra.setPosition( pos );
//		
//		RealTransformRandomAccess xfmRa = res.randomAccess();
//		int[] intPos = new int[]{10,10,10};
//		xfmRa.setPosition(intPos);
//		
//		
////		System.out.println(" xfm value at: (" + ArrayUtil.printArray(intPos) + "): " + val );
//		
////		System.out.println(" xfm value at: (" + ArrayUtil.printArray(intPos) + "): " + val );
//		
//		
////		double xfmval = rtra.get()ed
////		System.out.println(" xfm value at: (" + ArrayUtil.printArray(pos) + "): " + val );
//
//	}

	public static void testEdgelResamp2(){

		System.out.println("test Edgel Resamp");

		//Img<FloatType> testImg = ImgUtil.createGradientImgX( 20, 20, 20, new FloatType() );
		//Img<FloatType> testImg = ImgUtil.createGradientImgY( 20, 20, 20, new FloatType() );
		Img<FloatType> testImg = ImgUtil.createGradientImgZ( 20, 20, 20, new FloatType() );

		CrackCorrection<FloatType, UnsignedByteType> cc = 
				new CrackCorrection<FloatType, UnsignedByteType>( 
						testImg, null );
	
		//float[] norm = new float[]{ 1.0f, 1.0f, 0.0f };
		float[] norm = new float[]{ 0.0f, 0.0f, 1.0f };
		norm = ArrayUtil.normalizeLength(norm);

		System.out.println(" normal vector: (" + ArrayUtil.printArray(norm) + ") "  );

		Edgel edgel = new Edgel( new float[]{ 10.5f, 10.5f, 10.5f}, 
								 norm, 1);

		
		int[] patchSize = new int[]{5, 5, 3};
		Img<FloatType> resImg = cc.edgelToImage(edgel, testImg, patchSize);
		
		//ImageJFunctions.show(resImg);

	}
	
	public static void testEdgelResamp(){

		System.out.println("test Edgel Resamp");

		//Img<FloatType> testImg = ImgUtil.createGradientImgX( 20, 20, 20, new FloatType() );
		//Img<FloatType> testImg = ImgUtil.createGradientImgY( 20, 20, 20, new FloatType() );
		Img<FloatType> testImg = ImgUtil.createGradientImgZ( 20, 20, 20, new FloatType() );

		System.out.println(" testImg: " + testImg );
	
		//float[] norm = new float[]{ 1.0f, 1.0f, 0.0f };
		float[] norm = new float[]{ 0.0f, 0.0f, 1.0f };
		norm = ArrayUtil.normalizeLength(norm);

		System.out.println(" normal vector: (" + ArrayUtil.printArray(norm) + ") "  );

		Edgel edgel = new Edgel( new float[]{ 10.5f, 10.5f, 10.5f}, 
								 norm, 1);
		System.out.println("edgel : " + edgel);

		NLinearInterpolatorFactory<FloatType> interpFactory = 
			new NLinearInterpolatorFactory<FloatType>();

		RealRandomAccessible<FloatType> interpolant = 
			Views.interpolate( 
					Views.extendMirrorSingle( testImg ), interpFactory
			);

		RealRandomAccess<FloatType> rra = interpolant.realRandomAccess();		
		double[] pos = new double[]{ 10.5, 10.5, 10.5 } ;
		rra.setPosition( pos );
		double val = rra.get().getRealDouble();		

		System.out.println(" value at: (" + ArrayUtil.printArray(pos) + "): " + val );

		AffineTransform3D xfm = pickTransformation(edgel);
		
		int[] N = new int[]{5, 5, 3};
		Img<FloatType> resImg = normalPatch( pos, N, xfm, testImg );
		
		//ImageJFunctions.show(resImg);

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
	public static <T extends RealType<T>> Img<T> normalPatch(double[] basepos, int[] patchSize,  
			AffineTransform3D xfmIn,  Img<T> src ){
		
		AffineTransform3D xfm = xfmIn.copy();
		
		NLinearInterpolatorFactory< T > interpFactory =
		            new NLinearInterpolatorFactory< T >();
		 
		RealRandomAccessible< T > interpolant1 = Views.interpolate(
	            Views.extendMirrorSingle( src ), interpFactory );
		RealRandomAccess<T> srcRealRa = interpolant1.realRandomAccess();
		
		int ndims_in  = src.numDimensions();
		int ndims_out = patchSize.length;
		
		Img<T> out = src.factory().create( patchSize, Util.getTypeFromRandomAccess(src));
		
		logger.debug(" basepos  :" + ArrayUtil.printArray( basepos ) );
		logger.debug(" out size :" + ArrayUtil.printArray( patchSize ) );

		
		int[] 	 pos    = new int   [ndims_in];
		double[] dpos	= new double[ndims_in];
		
//		double[][][] a = new double[patchSize[0]][patchSize[1]][patchSize[2]];
		
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
		
//		double[][] c = new double[ndims_in][ndims_in+1];
//		xfm.toMatrix(c);
//		logger.debug("xfm: " + ArrayUtil.printArray(c));
		
		Cursor<T> cursor = out.cursor();
		while(cursor.hasNext()){
			
			cursor.fwd();
			cursor.localize(pos);
			
			double[] pt = ArrayUtil.toDouble(pos);
			
			xfm.apply( pt, dpos);
			
			logger.debug(" pt     :" + ArrayUtil.printArray( pt  ) );
			logger.debug(" dpos   :" + ArrayUtil.printArray( dpos) );
			
			
			srcRealRa.setPosition(dpos);
			double val = srcRealRa.get().getRealDouble();
			cursor.get().setReal( val );
			
			logger.debug(" val    :" + val );
			
//			a[pos[0]][pos[1]][pos[2]]  = val;
			
		}
		
//		System.out.println("res:\n" + ArrayUtil.printArray(a));
		
		return out;
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
		
		logger.debug(" Ra: \n" + ArrayUtil.printArray(Ra) );
		
		AffineTransform3D R = new AffineTransform3D();
		R.set(Ra);
		
		
//		AffineTransform3D R = RodriguesRotation.rotation(
//		new double[]{0,0,1}, ArrayUtil.toDouble(nrm));
		
		return R;
	}

   /**
	 * 
	 * @param in an M x N matrix whose rows span a subspace of R^N
	 * @return an orthogonal basis for the subspace of R^N not spanned by in
	 */
	public static DenseMatrix64F remainderSubspace(DenseMatrix64F in){
		SimpleSVD<SimpleMatrix> svd = new SimpleSVD<SimpleMatrix>(in,false);
		SimpleMatrix nullspace = svd.nullSpace();
		return nullspace.getMatrix();
	}


	   
	public static void main(String[] args){

//		testDistanceBasedCrackSideComp();	
		
//		testEdgelResamp();
		
//		testEdgelResamp2();
		
		testCrackCorrPipeline();	

		System.out.println("crack correction finished");
		System.exit(0);
	}


}
