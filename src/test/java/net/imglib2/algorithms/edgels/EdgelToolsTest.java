package net.imglib2.algorithms.edgels;

import static org.junit.Assert.*;


import net.imglib2.RandomAccess;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithms.edge.EdgelTools;
import net.imglib2.algorithms.patch.PatchTools;
import net.imglib2.img.Img;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.InverseRealTransform;
import net.imglib2.realtransform.RealTransformRandomAccessible;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgUtil;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.junit.Test;

import edu.jhu.ece.iacl.utility.ArrayUtil;

public class EdgelToolsTest {
	
	public static double tol = 0.01;
	
	public static Logger logger = LogManager.getLogger(EdgelToolsTest.class.getName());
	
	
	
	@Test
	public void testXfm() throws InterruptedException{
		
		Edgel e = new Edgel(
				new float[]{7f, 9f, 11f},
				new float[]{0f, 1f, 0f },
				1f
			);
			
		int[] patchSize = new int[]{5,5,3};
		int[] midPt = PatchTools.patchSizeToMidpt(patchSize);
		
		int[] expectMidPt = new int[]{2,2,1};
		
		assertArrayEquals("midpt 1", expectMidPt, midPt);
		
		AffineTransform3D xfm = EdgelTools.edgelToXfm(e, expectMidPt);
		
		logger.debug("xfm : " + xfm );
		
		double[] res = new double[3];
		double[] expectedResult000 = new double[]{7,9,11};
		xfm.apply( ArrayUtil.toDouble(expectMidPt), res);
		
		assertEquals(" xfm (0,0,0) x", expectedResult000[0], res[0], tol);
		assertEquals(" xfm (0,0,0) y", expectedResult000[1], res[1], tol );
		assertEquals(" xfm (0,0,0) z", expectedResult000[2], res[2], tol );

		double[] midPlusZ  = ArrayUtil.toDouble(ArrayUtil.clone(expectMidPt));
		midPlusZ[2]++;
		
		double[] midMinusZ = ArrayUtil.toDouble(ArrayUtil.clone(expectMidPt));
		midMinusZ[2]--;
		
		
		xfm.apply( midPlusZ, res);
		double[] expectMidPlusZ = new double[]{7,10,11};
		assertEquals(" xfm (0,0,1) x", expectMidPlusZ[0], res[0], tol);
		assertEquals(" xfm (0,0,1) y", expectMidPlusZ[1], res[1], tol );
		assertEquals(" xfm (0,0,1) z", expectMidPlusZ[2], res[2], tol );
		
		
		xfm.apply( midMinusZ, res);
		double[] expectMidMinusZ = new double[]{7,8,11};
		assertEquals(" xfm (0,0,-1) x", expectMidMinusZ[0], res[0], tol);
		assertEquals(" xfm (0,0,-1) y", expectMidMinusZ[1], res[1], tol );
		assertEquals(" xfm (0,0,-1) z", expectMidMinusZ[2], res[2], tol );
		
		
//		AffineTransform3D xfmInv = xfm.inverse();
		
	}
	
	@Test
	public void testView(){
		
		/**  test view **/
//		int nLevels = 5;
//		Img<FloatType> img = ImgUtil.createCheckerImg( new int[]{15,15,15}, new FloatType(), nLevels);
//		ImgUtil.writeFloat(img, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/checkerImg.tif");
		
		Img<FloatType> img = ImgUtil.createGradientImgY(15, 15, 15,  new FloatType());

		RandomAccess<FloatType> imgRa = img.randomAccess();
		
		int[] patchSize = new int[]{5,5,3};
		Edgel e = new Edgel(
				new float[]{7f, 9f, 11f},
				new float[]{0f, 1f, 0f },
				1f
			);
		
		RealTransformRandomAccessible<FloatType, InverseRealTransform> view = EdgelTools.edgelToView(e, img, patchSize);
		RealTransformRandomAccessible<FloatType, InverseRealTransform>.RealTransformRandomAccess viewRa = view.randomAccess();

		viewRa.setPosition(new int[]{2,2,0});
		imgRa.setPosition(new int[]{7,8,11});
		assertEquals("value at patch -z", imgRa.get().get(), viewRa.get().get(), tol);
		
		viewRa.setPosition(new int[]{2,2,1});
		imgRa.setPosition(new int[]{7,9,11});
		assertEquals("value at patch center", imgRa.get().get(), viewRa.get().get(), tol);
		
		viewRa.setPosition(new int[]{2,2,2});
		imgRa.setPosition(new int[]{7,10,11});
		assertEquals("value at patch +z", imgRa.get().get(), viewRa.get().get(), tol);
		
		
	}


}
