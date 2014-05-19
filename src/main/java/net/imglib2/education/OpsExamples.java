package net.imglib2.education;

import java.util.BitSet;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import ij.IJ;

import net.imglib2.histogram.Histogram1d;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.ops.function.real.RealImageFunction;
import net.imglib2.ops.function.real.StatCalculator;
import net.imglib2.ops.operation.iterableinterval.unary.MakeHistogram;
import net.imglib2.ops.pointset.IterableIntervalPointSet;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;

public class OpsExamples {

	static Logger logger = LogManager.getLogger( OpsExamples.class.getName() );

	public static <T extends RealType<T> & NativeType<T>> void imageHistogram( Img<T> img )
	{
		
		MakeHistogram<T> mh = new MakeHistogram<T>(256, true);
		Histogram1d<T> hist = mh.compute( img );

	}

	
	public static <T extends RealType<T> & NativeType<T>> void imageTextureStats( Img<T> img )
	{
		BitSet features2Comp = new BitSet(14);
		

	}

	public static <T extends RealType<T> & NativeType<T>> void imageStats( Img<T> img  )
	{

		RealImageFunction<T,T> imgfun = new RealImageFunction<T,T>( img, img.firstElement() );
		IterableIntervalPointSet ptset = new IterableIntervalPointSet( img );

		StatCalculator<T> sc = new	StatCalculator<T>( imgfun, ptset );

		double am = sc.arithmeticMean();
		double sd = sc.sampleStdDev();
		double sk = sc.sampleSkew();
		double kr = sc.sampleKurtosis();

		logger.info("am: " + am);
		logger.info("sd: " + sd);
		logger.info("sk: " + sk);
		logger.info("kr: " + kr);

	}

	public static void main( String[] args )
	{
		//String srcImgFn = "/groups/jain/home/bogovicj/learning/advanced-imglib2/images/bee-1.tif";
		//Img<FloatType> img =  ImagePlusAdapter.convertFloat( IJ.openImage(srcImgFn) );
	
		Img<FloatType> img = ImgOps.createGradientImgX(32,32,4, new FloatType());			

		imageStats( img );

		System.out.println("finished");
		System.exit(0);
	}
}
