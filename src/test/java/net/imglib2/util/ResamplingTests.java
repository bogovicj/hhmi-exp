package net.imglib2.util;

import ij.IJ;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class ResamplingTests {

	public static void main(String[] args) {
		
//		test1();
		test2();
		
	}

	public static void test1()
	{
		int[] sz = new int[]{ 21, 21, 21 };
		double[] dsFactors = new double[]{ 1, 1, 2 };
		
		Img<FloatType> img = ImgOps.createGradientImgX(sz, new FloatType());
		double[] sigmas = new double[]{ 0, 0, 0.5};
//		double[] sigmas = new double[]{ 2, 2, 0.5};
		
//		RandomAccessibleInterval<FloatType> view = Views.interval( 
//				Views.extendZero(img),
//				img);
		
		Img<FloatType> out = Resampling.resampleGaussian(
				img, 
				img.factory(), 
				dsFactors, sigmas, sigmas);
		
		System.out.println("out: " + out );
		
		int x = 11;
		int y = 11;
		
		RandomAccess<FloatType> ra = out.randomAccess();
		
		for( int z = 0; z < 11; z++ )
		{
			ra.setPosition(x, 0);
			ra.setPosition(y, 1);
			ra.setPosition(z, 2);
			
			System.out.println( " z " + z + "  " + ra.get());
		}
		
		ImgOps.writeFloat( out,  "/groups/saalfeld/home/bogovicj/tmp/grad_ds2.tif");
//		ImageJFunctions.show( out );
		
	}
	
	public static void test2(){
		String imfn = "/groups/saalfeld/home/bogovicj/tests/testdat/boats.tif";
		Img<FloatType> img =  ImagePlusAdapter.convertFloat( IJ.openImage(imfn) );
		double[] sigmas = new double[]{ 0, 0 };
		double[] dsFactors = new double[]{ 1, 1 };
		
		Img<FloatType> out = Resampling.resampleGaussian(
				img, 
				img.factory(), 
				dsFactors, sigmas, sigmas);
		
		ImgOps.writeFloat( out,  "/groups/saalfeld/home/bogovicj/tests/testdat/boats_myds_1-1.tif");
	}
}
