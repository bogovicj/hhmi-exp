package net.imglib2.algorithms.register;

import static org.junit.Assert.*;

import org.junit.Test;

import net.imglib2.Cursor;
import net.imglib2.ExtendedRandomAccessibleInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithms.edge.EdgelTools;
import net.imglib2.algorithms.patch.PatchTools;
import net.imglib2.algorithms.registration.TransformTools;
import net.imglib2.img.Img;
import net.imglib2.img.imageplus.ImagePlusImg;
import net.imglib2.img.imageplus.ImagePlusImgFactory;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.InverseRealTransform;
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
	public void testPcaReg(){
		
		System.out.println("testPcaReg");
		
		int[] sz = new int[]{ 27,27,27 };
		double[] ctr = ArrayUtil.toDouble( PatchTools.patchSizeToMidpt( sz ));
		double[] sigs = new double[]{0.1, 1, 10};
//		double mint = -1;
//		double maxt = 1.5;
		double mint = -1;
		double maxt = 9e20;
		
		double[] ctr2 = new double[]{ctr[0], sz[1]-5, ctr[2]}; 
		double[] sigs2 = new double[]{0.1, 0.1, 0.1};
		
		ImagePlusImgFactory<FloatType> factory = new ImagePlusImgFactory<FloatType>();
		
		Img<FloatType> img = ImgOps.createGaussianEllipseImg(factory, sz, ctr, sigs, mint, maxt, new FloatType());
		Img<FloatType> img2 = ImgOps.createGaussianEllipseImg(factory, sz, ctr2, sigs2, mint, maxt, new FloatType());
		
		RandomAccess<FloatType> imgRa = img.randomAccess();
		RandomAccess<FloatType> img2Ra = img2.randomAccess();
		
		Img<FloatType> img3 = img.factory().create(img, img.firstElement());
		Cursor<FloatType> c = img3.cursor();
		while(c.hasNext()){
			
			c.fwd();
			imgRa.setPosition(c);
			img2Ra.setPosition(c);
			
			img2Ra.get().mul(0.5);
			
			c.get().add(imgRa.get());
			c.get().add(img2Ra.get());
		}
		
//		ImgOps.writeFloat(img, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/ellipseImg.tif");
//		ImgOps.writeFloat(img2, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/ellipseImg2.tif");
		
		
//		double[] q = new double[4];
//		double[][] R = new double[3][3];
//		LinAlgHelpers.quaternionFromAngleAxis(new double[]{0, 0, 1}, Math.PI/8, q);
//		LinAlgHelpers.quaternionToR(q, R);

		AffineTransform3D xfm = TransformTools.rotationCentered(2, Math.PI/8, ctr);
		RealTransformRandomAccessible<FloatType, InverseRealTransform> img2Xfm = TransformTools.xfmToView(xfm, img3, new NLinearInterpolatorFactory<FloatType>());

		ImagePlusImg<FloatType, ?> img2xfmout = factory.create(img2, img2.firstElement());
		ImgOps.copyInto( img2Xfm, img2xfmout);
		ImgOps.writeFloat(img2xfmout, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/ellipseImg2xfm.tif");
		System.out.println("done");
		
		assertTrue(true);
		
	}

}
