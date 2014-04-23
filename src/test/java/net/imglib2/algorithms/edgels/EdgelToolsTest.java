package net.imglib2.algorithms.edgels;

import static org.junit.Assert.*;

import java.util.ArrayList;

import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithms.crack.CrackCorrection;
import net.imglib2.algorithms.edge.EdgelTools;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.InverseRealTransform;
import net.imglib2.realtransform.RealTransformRandomAccessible;
import net.imglib2.realtransform.RealTransformRandomAccessible.RealTransformRandomAccess;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgUtil;
import net.imglib2.view.Views;

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
		int[] midPt = EdgelTools.patchSizeToMidpt(patchSize);
		
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
		
		// test view
		Img<FloatType> img = ImgUtil.createGradientImgY(32, 32, 32, new FloatType());
		RealTransformRandomAccessible<FloatType, InverseRealTransform> view = EdgelTools.edgelToView(e, img, patchSize);
		RealTransformRandomAccess viewRa = view.randomAccess();
		
		viewRa.setPosition(new int[]{0,0,0});
		logger.info("val: " + viewRa.get());
		
		viewRa.setPosition(new int[]{0,0,-1});
		logger.info("val: " + viewRa.get());

		viewRa.setPosition(new int[]{0,0,1});
		logger.info("val: " + viewRa.get());
		
		ArrayImgFactory<UnsignedByteType> ubfactory = new ArrayImgFactory<UnsignedByteType>();
		Img<UnsignedByteType> maskImg = ubfactory.create(img, new UnsignedByteType());
		
		UnsignedByteType maskVal = new UnsignedByteType(1);
		
		CrackCorrection.setMask( e.getPosition(), 
				patchSize, xfm, Views.extendValue(maskImg, new UnsignedByteType(0)), maskVal);
		

		assertEquals( "window size ", patchSize[0]*patchSize[1]*patchSize[2], ImgUtil.numNonZero(maskImg));
		
		
	}
	

}
