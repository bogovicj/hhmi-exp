package net.imglib2.algorithms.edgels;

import static org.junit.Assert.*;

import org.ejml.data.DenseMatrix64F;
import org.junit.Test;

import java.util.Arrays;

import net.imglib2.RandomAccess;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithms.crack.CrackCorrection;
import net.imglib2.algorithms.edge.EdgelTools;
import net.imglib2.algorithms.patch.PatchTools;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.InverseRealTransform;
import net.imglib2.realtransform.RealTransformRandomAccessible;
import net.imglib2.realtransform.RealTransformRandomAccessible.RealTransformRandomAccess;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;
import edu.jhu.ece.iacl.utility.ArrayUtil;

public class EdgelSamplingTest {

	@Test
	public void testEdgelView(){

		// generate image with gradient in X
		Img<FloatType> img = ImgOps.createGradientImgX( 
				new int[]{32,32,32}, new FloatType());

		// choose an edgel with normal in X direction
		Edgel e = new Edgel( new double[]{12,12,12},
									new double[]{1,0,0}, 1);
		
		// expect a constant value in the edgel view
		int[] patchSize  = new int[]{5,5,3};

		//RealTransformRandomAccessible<FloatType,InverseRealTransform> view 
		RealTransformRandomAccessible<FloatType, InverseRealTransform> view =
				EdgelTools.edgelToView(e, img, patchSize);
		
		
		RealTransformRandomAccess vra = view.randomAccess();
		
		vra.setPosition(new int[]{0,0,0});	
		System.out.println("val (0,0,0) : " + vra.get());
		
		vra.setPosition(new int[]{0,1,0});	
		System.out.println("val (0,1,0) : " + vra.get());
		
		vra.setPosition(new int[]{1,0,0});	
		System.out.println("val (1,0,0) : " + vra.get());
		
		vra.setPosition(new int[]{1,1,0});	
		System.out.println("val (1,1,0) : " + vra.get());
		
		int[] midPt = PatchTools.patchSizeToMidpt(patchSize);
		AffineTransform3D xfmIn = EdgelTools.edgelToXfm(e, midPt);
		
		
		ArrayImgFactory<UnsignedByteType> ubfactory = new ArrayImgFactory<UnsignedByteType>();
		
		Img<UnsignedByteType> maskimg = ubfactory.create(img, new UnsignedByteType());
		
		double[] pos = new double[e.numDimensions()];
		e.localize(pos);
		CrackCorrection.setMask( pos, patchSize, xfmIn, maskimg, 
								 new UnsignedByteType(255));
		
//		ImgUtil.write(maskimg, "/groups/jain/home/bogovicj/projects/crackSegmentation/edgelValidate/toy.tif");
		
//		assertEquals("hi",1,1);
	}

}
