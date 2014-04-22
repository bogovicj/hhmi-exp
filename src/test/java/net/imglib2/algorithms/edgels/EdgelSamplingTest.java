package net.imglib2.algorithms.edgels;

import static org.junit.Assert.*;

import org.ejml.data.DenseMatrix64F;
import org.junit.Test;

import java.util.Arrays;

import net.imglib2.RandomAccess;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithms.crack.CrackCorrection;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.InverseRealTransform;
import net.imglib2.realtransform.RealTransformRandomAccessible;
import net.imglib2.realtransform.RealTransformRandomAccessible.RealTransformRandomAccess;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgUtil;
import edu.jhu.ece.iacl.utility.ArrayUtil;

public class EdgelSamplingTest {

	@Test
	public void testEdgelView(){

		// generate image with gradient in X
		Img<FloatType> img = ImgUtil.createGradientImgX(32,32,32, new FloatType(0));

		// choose an edgel with normal in X direction
		Edgel e = new Edgel( new float[]{12,12,12},
									new float[]{1f,0f,0f}, 1f);
		
		// expect a constant value in the edgel view
		int[] patchSize  = new int[]{5,5,3};

		//RealTransformRandomAccessible<FloatType,InverseRealTransform> view 
		RealTransformRandomAccessible<FloatType, InverseRealTransform> view =
				CrackCorrection.edgelToView(e, img, patchSize);
		
		RealTransformRandomAccess vra = view.randomAccess();
		
		vra.setPosition(new int[]{0,0,0});	
		System.out.println("val (0,0,0) : " + vra.get());
		
		vra.setPosition(new int[]{0,1,0});	
		System.out.println("val (0,1,0) : " + vra.get());
		
		vra.setPosition(new int[]{1,0,0});	
		System.out.println("val (1,0,0) : " + vra.get());
		
		vra.setPosition(new int[]{1,1,0});	
		System.out.println("val (1,1,0) : " + vra.get());
		
		int[] midPt = CrackCorrection.patchSizeToMidpt(patchSize);
		AffineTransform3D xfmIn = CrackCorrection.edgelToXfm(e, midPt);
		
		
		ArrayImgFactory<UnsignedByteType> ubfactory = new ArrayImgFactory<UnsignedByteType>();
		
		Img<UnsignedByteType> maskimg = ubfactory.create(img, new UnsignedByteType());
		
		CrackCorrection.setMask( ArrayUtil.toDouble(e.getPosition()), 
				patchSize, xfmIn, maskimg, new UnsignedByteType(255));
		
		ImgUtil.write(maskimg, "/groups/jain/home/bogovicj/projects/crackSegmentation/edgelValidate/toy.tif");
		
		
//		assertEquals("hi",1,1);
	}

}
