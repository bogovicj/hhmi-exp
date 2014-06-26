package net.imglib2.util;

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class ResamplingTests {

	public static void main(String[] args) {
		
		test1();
		
	}

	public static void test1()
	{
		int[] sz = new int[]{ 21, 21, 21 };
		double[] dsFactors = new double[]{ 1, 1, 2 };
		
		Img<FloatType> img = ImgOps.createGradientImgX(sz, new FloatType());
		double[] sigmas = new double[]{ 1, 1, 2};
		
//		RandomAccessibleInterval<FloatType> view = Views.interval( 
//				Views.extendZero(img),
//				img);
		
		Img<FloatType> out = Resampling.resampleGaussian(
				img, 
				img.factory(), 
				dsFactors, sigmas, sigmas);
		
		System.out.println("out: " + out );
		
	}
}
