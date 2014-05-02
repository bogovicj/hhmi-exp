package net.imglib2.algorithms.register;

import static org.junit.Assert.*;

import org.apache.log4j.Logger;
import org.ejml.data.DenseMatrix64F;
import org.junit.Test;

import net.imglib2.Cursor;
import net.imglib2.ExtendedRandomAccessibleInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithms.edge.EdgelTools;
import net.imglib2.algorithms.moments.ImgMoment;
import net.imglib2.algorithms.patch.PatchTools;
import net.imglib2.algorithms.registration.TransformTools;
import net.imglib2.img.Img;
import net.imglib2.img.imageplus.ImagePlusImg;
import net.imglib2.img.imageplus.ImagePlusImgFactory;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.realtransform.AffineTransform;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.InverseRealTransform;
import net.imglib2.realtransform.RealTransform;
import net.imglib2.realtransform.RealTransformRandomAccessible;
import net.imglib2.realtransform.RealTransformSequence;
import net.imglib2.realtransform.RealViews;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;
import net.imglib2.util.LinAlgHelpers;
import net.imglib2.view.Views;
import edu.jhu.ece.iacl.utility.ArrayUtil;

public class TestPcaRotReg {
	

	@Test
	public void testPca(){
	System.out.println("testPcaReg");
		
		int[] sz = new int[]{ 27,27,27 };
		double[] ctr = ArrayUtil.toDouble( PatchTools.patchSizeToMidpt( sz ));
		double[] sigs = new double[]{1, 10, 0.1};
//		double mint = -1;
//		double maxt = 1.5;
		double mint = -1;
		double maxt = 9e20;
		
		ImagePlusImgFactory<FloatType> factory = new ImagePlusImgFactory<FloatType>();
		
		Img<FloatType> img = ImgOps.createGaussianEllipseImg(factory, sz, ctr, sigs, mint, maxt, new FloatType());
		ImgOps.writeFloat(img, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/ellipseImg.tif");

		AffineTransform3D xfm = TransformTools.rotationCentered(2, Math.PI/8, ctr);
		RealTransformRandomAccessible<FloatType, InverseRealTransform> imgXfm = 
				TransformTools.xfmToView(xfm, 
						Views.extendZero(img), new NLinearInterpolatorFactory<FloatType>());

		ImagePlusImg<FloatType, ?> imgXfmOut = factory.create(img, img.firstElement());
		ImgOps.copyInto( imgXfm, imgXfmOut);

		ImgMoment tgtMom = new ImgMoment();
		double[] tgtOr = tgtMom.orientation(imgXfmOut);
		tgtMom.orientationEvalsEvecs(tgtOr);
		
		System.out.println("tgt or: " + ArrayUtil.printArray(tgtOr));
		System.out.println("tgt cent: " + ArrayUtil.printArray(tgtMom.centroid()));
		
		
		DenseMatrix64F srcR = new DenseMatrix64F( 3, 3 );
		srcR.setData( tgtMom.getEvecs() );
		DenseMatrix64F srcMtx = TransformTools.rotationCenteredMtx(srcR, tgtMom.centroid());
	
	}
	
	
	@Test
	public void testPcaReg(){
		
		System.out.println("testPcaReg");
		
		int[] sz = new int[]{ 27,27,27 };
		double[] ctr = ArrayUtil.toDouble( PatchTools.patchSizeToMidpt( sz ));
		double[] sigs = new double[]{1, 10, 0.1};
//		double mint = -1;
//		double maxt = 1.5;
		double mint = -1;
		double maxt = 9e20;
		
		ImagePlusImgFactory<FloatType> factory = new ImagePlusImgFactory<FloatType>();
		
		Img<FloatType> img = ImgOps.createGaussianEllipseImg(factory, sz, ctr, sigs, mint, maxt, new FloatType());
		ImgOps.writeFloat(img, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/ellipseImg.tif");

		AffineTransform3D xfm = TransformTools.rotationCentered(2, Math.PI/8, ctr);
		RealTransformRandomAccessible<FloatType, InverseRealTransform> imgXfm = 
				TransformTools.xfmToView(xfm, 
						Views.extendZero(img), new NLinearInterpolatorFactory<FloatType>());

		ImagePlusImg<FloatType, ?> imgXfmOut = factory.create(img, img.firstElement());
		ImgOps.copyInto( imgXfm, imgXfmOut);
		ImgOps.writeFloat(imgXfmOut, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/ellipseImgXfm.tif");
		
		
		AffineTransform xfmEst = TransformTools.rotationPca( imgXfmOut, img);
		System.out.println("xfmEst:\n" + TransformTools.printAffineTransform(xfmEst));
		
//		AffineTransform xfmCat = xfmEst.preConcatenate(xfm);
//		System.out.println("xfmCat:\n" + TransformTools.printAffineTransform(xfmCat));
		
		
		RealTransformRandomAccessible<FloatType, InverseRealTransform> img2XfmXfm = 
				TransformTools.xfmToView( xfmEst, 
						Views.extendZero(imgXfmOut), new NLinearInterpolatorFactory<FloatType>());
		
		ImagePlusImg<FloatType, ?> img2xfmxfmout = factory.create(img, img.firstElement());
		ImgOps.copyInto( img2XfmXfm, img2xfmxfmout);
		ImgOps.writeFloat(img2xfmxfmout, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/ellipseImg2XfmInv.tif");
		
		
		System.out.println("done");
		assertTrue(true);
		
	}
	
	@Test
	public void testPcaReg2(){
		
		System.out.println("testPcaReg");
		
		int[] sz = new int[]{ 27,27,27 };
		double[] ctr = ArrayUtil.toDouble( PatchTools.patchSizeToMidpt( sz ));
		double[] sigs = new double[]{10, 1, 0.1};
//		double mint = -1;
//		double maxt = 1.5;
		double mint = -1;
		double maxt = 9e20;
		
		double[] ctrx = new double[]{sz[0]-5, ctr[1], ctr[2]};
		double[] ctry = new double[]{ctr[0], sz[1]-5, ctr[2]};
		double[] ctrz = new double[]{ctr[0], ctr[1], sz[2]-5};
		double[] sigs2 = new double[]{0.1, 0.1, 0.1};
		
		ImagePlusImgFactory<FloatType> factory = new ImagePlusImgFactory<FloatType>();
		
		Img<FloatType> img = ImgOps.createGaussianEllipseImg(factory, sz, ctr, sigs, mint, maxt, new FloatType());
		
		Img<FloatType> bx = ImgOps.createGaussianEllipseImg(factory, sz, ctrx, sigs2, mint, maxt, new FloatType());
		Img<FloatType> by = ImgOps.createGaussianEllipseImg(factory, sz, ctry, sigs2, mint, maxt, new FloatType());
		Img<FloatType> bz = ImgOps.createGaussianEllipseImg(factory, sz, ctrz, sigs2, mint, maxt, new FloatType());
		
		RandomAccess<FloatType> imgRa = img.randomAccess();
		
		RandomAccess<FloatType> bxRa = bx.randomAccess();
		RandomAccess<FloatType> byRa = by.randomAccess();
		RandomAccess<FloatType> bzRa = bz.randomAccess();
		
		Img<FloatType> img3 = img.factory().create(img, img.firstElement());
		Cursor<FloatType> c = img3.cursor();
		while(c.hasNext()){
			
			c.fwd();
			imgRa.setPosition(c);
			bxRa.setPosition(c);
			byRa.setPosition(c);
			bzRa.setPosition(c);
			
			bxRa.get().mul(0.2);
			byRa.get().mul(0.1);
			bzRa.get().mul(0.05);
			
			c.get().add(imgRa.get());
			c.get().add(bxRa.get());
			c.get().add(byRa.get());
			c.get().add(bzRa.get());
		}
		
		ImgOps.writeFloat(img3, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/ellipseImg3.tif");
		
//		ImgOps.writeFloat(img, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/ellipseImg.tif");
//		ImgOps.writeFloat(img2, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/ellipseImg2.tif");
		
//		double[] q = new double[4];
//		double[][] R = new double[3][3];
//		LinAlgHelpers.quaternionFromAngleAxis(new double[]{0, 0, 1}, Math.PI/8, q);
//		LinAlgHelpers.quaternionToR(q, R);

		AffineTransform3D xfm = TransformTools.rotationCentered(2, Math.PI/8, ctr);
		RealTransformRandomAccessible<FloatType, InverseRealTransform> img3Xfm = 
				TransformTools.xfmToView(xfm, Views.extendZero(img3), new NLinearInterpolatorFactory<FloatType>());

		ImagePlusImg<FloatType, ?> img3xfmout = factory.create(img, img.firstElement());
		ImgOps.copyInto( img3Xfm, img3xfmout);
		ImgOps.writeFloat(img3xfmout, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/ellipseImg3Xfm.tif");
		
		
		AffineTransform xfmEst = TransformTools.rotationPca( img3xfmout, img);
		System.out.println("xfmEst:\n" + TransformTools.printAffineTransform(xfmEst));
		
//		AffineTransform xfmCat = xfmEst.preConcatenate(xfm);
//		System.out.println("xfmCat:\n" + TransformTools.printAffineTransform(xfmCat));
		
		RealTransformRandomAccessible<FloatType, InverseRealTransform> img3XfmXfm = 
				TransformTools.xfmToView(xfmEst, Views.extendZero(img3xfmout), new NLinearInterpolatorFactory<FloatType>());
		
		ImagePlusImg<FloatType, ?> img3xfmxfmout = factory.create(img, img.firstElement());
		ImgOps.copyInto( img3XfmXfm, img3xfmxfmout);
		ImgOps.writeFloat(img3xfmxfmout, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/ellipseImg3XfmInv.tif");
		
		
		System.out.println("done");
		assertTrue(true);
		
	}

}
