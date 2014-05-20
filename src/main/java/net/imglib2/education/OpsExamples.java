package net.imglib2.education;

import java.util.BitSet;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import ij.IJ;
import ij.ImageJ;
import net.imglib2.histogram.Histogram1d;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.ops.condition.FunctionNotEqualCondition;
import net.imglib2.ops.data.CooccurrenceMatrix;
import net.imglib2.ops.function.Function;
import net.imglib2.ops.function.real.RealImageFunction;
import net.imglib2.ops.function.real.RealSumFunction;
import net.imglib2.ops.function.real.StatCalculator;
import net.imglib2.ops.img.SerialImageAssignment;
import net.imglib2.ops.img.UnaryOperationAssignment;
import net.imglib2.ops.input.PointInputIterator;
import net.imglib2.ops.input.PointSetInputIterator;
import net.imglib2.ops.operation.UnaryOperation;
import net.imglib2.ops.operation.img.binary.ImgCombine;
import net.imglib2.ops.operation.iterableinterval.unary.MakeCooccurrenceMatrix;
import net.imglib2.ops.operation.iterableinterval.unary.MakeCooccurrenceMatrix.HaralickFeature;
import net.imglib2.ops.operation.iterableinterval.unary.MakeHistogram;
import net.imglib2.ops.operation.real.binary.RealAdd;
import net.imglib2.ops.operation.real.binary.RealBinaryOperation;
import net.imglib2.ops.operation.real.unary.Normalize;
import net.imglib2.ops.pointset.IterableIntervalPointSet;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;
import net.imglib2.view.Views;

public class OpsExamples {

	static Logger logger = LogManager.getLogger(OpsExamples.class.getName());

	public static <T extends RealType<T> & NativeType<T>> void imageHistogram(
			Img<T> img) {

		MakeHistogram<T> mh = new MakeHistogram<T>(256, true);
		Histogram1d<T> hist = mh.compute(img);

	}

	public static <T extends RealType<T> & NativeType<T>> void imageTextureStats(
			Img<T> img) {
		BitSet features2Comp = new BitSet(14);
		features2Comp.set(0, 14);

		int nValues = 32;

		CooccurrenceMatrix coMtx = new CooccurrenceMatrix(nValues);
		MakeCooccurrenceMatrix<T> mc = new MakeCooccurrenceMatrix<T>(0, 1, 1,
				nValues, CooccurrenceMatrix.MatrixOrientation.HORIZONTAL,
				features2Comp);

		mc.compute(img, coMtx);

		logger.info("coMtx: " + coMtx);

		HaralickFeature[] values = MakeCooccurrenceMatrix.HaralickFeature.values();
		for( int i=0; i<values.length; i++)
		{
			logger.info("feature " + values[i] + " value: " + coMtx.getFeature(i) );
		}

	}

	public static <T extends RealType<T> & NativeType<T>> void evalFun( Img<T> img1, Img<T> img2, RealBinaryOperation<T,T,T> fun  )
	{
		Img<T> res = img1.factory().create( img1, img1.firstElement() );
		RealImageFunction<T,T> imgfun = new RealImageFunction<T,T>( img1, img1.firstElement() );
		
		ImgCombine<T,T,T> imcomb = new ImgCombine<T,T,T>(fun);
		imcomb.compute(img1, img2, res);
		
		ImgOps.writeFloat(res, "/groups/jain/home/bogovicj/learning/advanced-imglib2/images/bee-1-wave.tif");
		
	}
	
	public static <T extends RealType<T> & NativeType<T>> void evalFun( Img<T> img, UnaryOperation<T,T> fun  )
	{
		Img<T> res = img.factory().create( img, img.firstElement() );

		
		UnaryOperationAssignment<T,T> uoa = new UnaryOperationAssignment<T,T>(  fun) ;
		uoa.compute( img, res );
		
		ImgOps.writeFloat(res, "/groups/jain/home/bogovicj/learning/advanced-imglib2/images/bee-1-nrm.tif");
		
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
		String srcImgFn = "/groups/jain/home/bogovicj/learning/advanced-imglib2/images/bee-1.tif";
		Img<FloatType> img =  ImagePlusAdapter.convertFloat( IJ.openImage(srcImgFn) );
	
		int[] sz = new int[]{ (int)img.dimension(0), (int)img.dimension(1) };
		Img<FloatType> img2 = ImgOps.createGradientImgX( sz, new FloatType());			

		//imageStats( img );
		//imageTextureStats( img );
	
		// function that normalizes to range [0,1] 
//		Normalize<FloatType> nrmFun = new Normalize<FloatType>( 
//			0, 255, 0, 1);	
//		evalFun( img, nrmFun );

		// add the images
		RealAdd<FloatType,FloatType,FloatType> sumfun = new RealAdd<FloatType,FloatType,FloatType>();
		evalFun( img, img2, sumfun );

		System.out.println("finished");
		System.exit(0);
	}
}
