package net.imglib2.algorithms.crack.exps;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import edu.jhu.ece.iacl.utility.ArrayUtil;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithm.edge.SubpixelEdgelDetection;
import net.imglib2.algorithms.crack.CrackCorrection;
import net.imglib2.exception.ImgLibException;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.imageplus.ImagePlusImg;
import net.imglib2.img.imageplus.ImagePlusImgs;
import net.imglib2.io.ImgIOException;
import net.imglib2.io.ImgOpener;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgUtil;
import net.imglib2.util.LinAlgHelpers;
import net.imglib2.view.Views;

public class EdgelOrientationVisValid {

	protected static Logger logger = LogManager.getLogger(EdgelOrientationVisValid.class
			.getName());
	
	public static void visualizeEdgelBoxes() 
	{
		
		String imgfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/crackVolDown_cp.tif";
		
//		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp.tif";
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_smooth.tif";

		String edgelMaskPrefix = "/groups/jain/home/bogovicj/projects/crackSegmentation/edgelValidate/edgelMask_1";

		ArrayImgFactory<UnsignedByteType> ubfactory = new ArrayImgFactory<UnsignedByteType>();
		ArrayImgFactory<IntType> ifactory = new ArrayImgFactory<IntType>();
		ArrayImgFactory<FloatType> ffactory = new ArrayImgFactory<FloatType>();
		
		int[] patchSize = new int[]{7,7,3};
		int[] patchMidPt = CrackCorrection.patchSizeToMidpt(patchSize);
		

		Img<FloatType> img = null;
//		Img<UnsignedByteType> mask = null;
		Img<FloatType> mask = null;
		try {
			img = new ImgOpener().openImg(imgfn, ffactory, new FloatType());
			
//			mask = new ImgOpener().openImg(maskfn, ubfactory,
//					new UnsignedByteType());
			mask = new ImgOpener().openImg(maskfn, ffactory,
					new FloatType());
			
		} catch (ImgIOException e) {
			e.printStackTrace();
		}

//		CrackCorrection<FloatType, UnsignedByteType> cc = new CrackCorrection<FloatType, UnsignedByteType>(
//				img, mask);
		CrackCorrection<FloatType, FloatType> cc = new CrackCorrection<FloatType, FloatType>(
				img, mask);
		
		cc.computeEdgels();
		
		int N = 40;
		
		ArrayList<Edgel> edgels = cc.getEdgels();
		Collections.shuffle(edgels, new Random(1));
		
		Img<UnsignedByteType> edgelmaskimg =ubfactory.create(mask, new UnsignedByteType(0));
		
		int i = 0;
		
		while(i < N){
			
			
			Edgel edgel = edgels.get(i);
			
//			IntType val = new IntType( i + 1 );
			UnsignedByteType val = new UnsignedByteType( 255 - i  );
			
//			AffineTransform3D xfmIn = CrackCorrection.pickTransformation(edgel);
//			CrackCorrection.edgelToView(edgel, edgelmaskimg, patchSize);
			
			
			AffineTransform3D xfmIn = CrackCorrection.edgelToXfm(edgel, patchMidPt);
			
			CrackCorrection.setMask( ArrayUtil.toDouble(edgel.getPosition()), 
					patchSize, xfmIn, 
					Views.extendValue(edgelmaskimg, new UnsignedByteType()), val);

			i++;
		}
		
		ImgUtil.write(edgelmaskimg, edgelMaskPrefix + ".tif");
		
		float[] pos = new float[]{221,110,6};
		int j = cc.edgelIdxNearest(pos);
		
		Edgel edgel = edgels.get(j);
		
		logger.info(" edgel pos : " + ArrayUtil.printArray(edgel.getPosition()));
		logger.info(" edgel grad: " + ArrayUtil.printArray(edgel.getGradient()));
		logger.info(" edgel mag : " + edgel.getMagnitude());
		
//		int[] coord = new int[3];
//		
//		RandomAccess<UnsignedByteType> ra = mask.randomAccess();
//		
//		for(int x=-1; x<1; x++)for(int y=-1; y<1; y++)for(int z=-1; z<1; z++){
//			
//			coord[0] = (int)Math.round( edgel.getPosition()[0] + x);
//			coord[1] = (int)Math.round( edgel.getPosition()[1] + y);
//			coord[2] = (int)Math.round( edgel.getPosition()[2] + z);
//			
//			logger.info(" pos " + ArrayUtil.printArray( coord ));
//			
//			
//			
//		}
		
	}
	
	public static void checkEdgelPosition() 
	{
		
		String imgfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/crackVolDown_cp.tif";
		
//		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp.tif";
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_smooth.tif";

		String edgelMaskPrefix = "/groups/jain/home/bogovicj/projects/crackSegmentation/edgelValidate/edgelMask_1";

		ArrayImgFactory<UnsignedByteType> ubfactory = new ArrayImgFactory<UnsignedByteType>();
		ArrayImgFactory<IntType> ifactory = new ArrayImgFactory<IntType>();
		ArrayImgFactory<FloatType> ffactory = new ArrayImgFactory<FloatType>();
		
		int[] patchSize = new int[]{7,7,3};
		int[] patchMidPt = CrackCorrection.patchSizeToMidpt(patchSize);
		
		Img<FloatType> img = null;
//		Img<UnsignedByteType> mask = null;
		Img<FloatType> mask = null;
		try {
			img = new ImgOpener().openImg(imgfn, ffactory, new FloatType());
			
//			mask = new ImgOpener().openImg(maskfn, ubfactory,
//					new UnsignedByteType());
			mask = new ImgOpener().openImg(maskfn, ffactory,
					new FloatType());
			
		} catch (ImgIOException e) {
			e.printStackTrace();
		}

		CrackCorrection<FloatType, FloatType> cc = new CrackCorrection<FloatType, FloatType>(
				img, mask);
		
		cc.computeEdgels();
		ArrayList<Edgel> edgels = cc.getEdgels();
		
		float[] pos = new float[]{206,137,9};
		int i = cc.edgelIdxNearest(pos);
		
		Edgel edgel = edgels.get(i);
		
		logger.info(" edgel pos : " + ArrayUtil.printArray(edgel.getPosition()));
		logger.info(" edgel grad: " + ArrayUtil.printArray(edgel.getGradient()));
		logger.info(" edgel mag : " + edgel.getMagnitude());
		
		
	}
	
	public static void checkEdgelType() 
	{
		
		String imgfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/crackVolDown_cp.tif";
		
//		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp.tif";
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_smooth.tif";

		String edgelMaskPrefix = "/groups/jain/home/bogovicj/projects/crackSegmentation/edgelValidate/edgelMask_1";
		String crackMaskRewriteFn = "/groups/jain/home/bogovicj/projects/crackSegmentation/edgelValidate/crackMaskOut.tif";

		ArrayImgFactory<UnsignedByteType> ubfactory = new ArrayImgFactory<UnsignedByteType>();
		ArrayImgFactory<IntType> ifactory = new ArrayImgFactory<IntType>();
		ArrayImgFactory<FloatType> ffactory = new ArrayImgFactory<FloatType>();
		
		int[] patchSize = new int[]{7,7,3};
		int[] patchMidPt = CrackCorrection.patchSizeToMidpt(patchSize);
		
		Img<FloatType> mask = null;
//		Img<UnsignedByteType> mask = null;
		try {
			
			mask = new ImgOpener().openImg(maskfn, ffactory,
					new FloatType());
			
//			mask = new ImgOpener().openImg(maskfn, ubfactory,
//					new UnsignedByteType());
			
			
		} catch (ImgIOException e) {
			e.printStackTrace();
		}

		ArrayList<Edgel> edgels = SubpixelEdgelDetection.getEdgels(mask, mask.factory(), 10.0f);
		logger.info("num edgels " + edgels.size());
		
		
		float[] pos = new float[]{206,137,9};

		
//		ImgUtil.write(mask, crackMaskRewriteFn);
		CrackCorrection cc = new CrackCorrection();
		cc.setEdgels(edgels);
		
		int i = cc.edgelIdxNearest(pos);
		
		Edgel edgel = edgels.get(i);
//		
		logger.info(" edgel pos : " + ArrayUtil.printArray(edgel.getPosition()));
		logger.info(" edgel grad: " + ArrayUtil.printArray(edgel.getGradient()));
		logger.info(" edgel mag : " + edgel.getMagnitude());
//		
		
	}
	
	public static void ball() throws ImgLibException{
		String edgelMaskPrefix = "/groups/jain/home/bogovicj/projects/crackSegmentation/edgelValidate/edgelball_1";
		String edgelPrefix = "/groups/jain/home/bogovicj/projects/crackSegmentation/edgelValidate/ball_1";
		
		ImagePlusImg<FloatType,?> img = ImagePlusImgs.floats(new long[]{200, 200, 200});
		Cursor<FloatType> c = Views.iterable(Views.offset(img, new long[]{100, 100, 100})).localizingCursor();
		double[] x = new double[]{0, 0, 0};
		while (c.hasNext()) {
		  FloatType t = c.next();
		  c.localize(x);
		  if (LinAlgHelpers.squareLength(x) < 2500)
		         t.set(0xff);
		}
//		img.getImagePlus().show();
		
		
		CrackCorrection<FloatType, FloatType> cc = new CrackCorrection<FloatType, FloatType>(
				img, img);
		
		cc.computeEdgels();
		
		int N = 40;
		
		ArrayList<Edgel> edgels = cc.getEdgels();
		Collections.shuffle(edgels, new Random(1));
		
		ArrayImgFactory<UnsignedByteType> ubfactory = new ArrayImgFactory<UnsignedByteType>();
		Img<UnsignedByteType> edgelmaskimg = ubfactory.create(img, new UnsignedByteType(0));
		
		int i = 0;
		
		int[] patchSize = new int[]{7,7,3};
		int[] patchMidPt = CrackCorrection.patchSizeToMidpt(patchSize);
		
		while(i < N){
			
			
			Edgel edgel = edgels.get(i);
			
//			IntType val = new IntType( i + 1 );
			UnsignedByteType val = new UnsignedByteType( 255 - i  );
			
//			AffineTransform3D xfmIn = CrackCorrection.pickTransformation(edgel);
//			CrackCorrection.edgelToView(edgel, edgelmaskimg, patchSize);
			
			
			AffineTransform3D xfmIn = CrackCorrection.edgelToXfm(edgel, patchMidPt);
			
			CrackCorrection.setMask( ArrayUtil.toDouble(edgel.getPosition()), 
					patchSize, xfmIn, 
					Views.extendValue(edgelmaskimg, new UnsignedByteType()), val);

			i++;
		}
		
		
		
		ImgUtil.write(edgelmaskimg, edgelMaskPrefix + ".tif");
		ImgUtil.write(img, edgelPrefix + ".tif");

		
//		ImgUtil.write(mask, crackMaskRewriteFn);
		
	}
	
	public static void checkImport(){
		
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_smooth.tif";
		
		ArrayImgFactory<FloatType> ffactory = new ArrayImgFactory<FloatType>();
		ArrayImgFactory<UnsignedByteType> ubfactory = new ArrayImgFactory<UnsignedByteType>();
		
		Img<FloatType> mask = null;
//		Img<UnsignedByteType> mask = null;
		try {
			
			mask = new ImgOpener().openImg(maskfn, ffactory,
					new FloatType());
			
//			mask = new ImgOpener().openImg(maskfn, ubfactory,
//					new UnsignedByteType());
			
			
		} catch (ImgIOException e) {
			e.printStackTrace();
		}
		
//		ArrayList<Integer> set = ImgUtil.uniqueInt(mask);
//		Collections.sort(set);
//		System.out.println("set: " + set);
		
		
//		RandomAccess<UnsignedByteType> ra = mask.randomAccess();
		RandomAccess<FloatType> ra = mask.randomAccess();
		
		ra.setPosition(new int[]{207,137,9});
		
		System.out.println("val: " + ra.get());
		
	}
	
	public static void main(String[] args) throws ImgLibException {
		System.out.println("starting");

		
//		visualizeEdgelBoxes();
//		checkEdgelPosition();
		
//		checkEdgelType();
//		checkImport();
		
		ball();
		
		System.out.println("finished");
//		System.exit(0);
	}

}
