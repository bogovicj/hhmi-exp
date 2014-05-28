package net.imglib2.algorithms.edgels;

import static org.junit.Assert.*;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithm.gradient.PartialDerivative;
import net.imglib2.algorithms.edge.EdgelTools;
import net.imglib2.algorithms.patch.PatchTools;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.iterator.IntervalIterator;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.InverseRealTransform;
import net.imglib2.realtransform.RealTransformRandomAccessible;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.junit.Test;

import edu.jhu.ece.iacl.utility.ArrayUtil;

public class EdgelToolsTest {
	
	public static double tol = 0.01;
	
	public static Logger logger = LogManager.getLogger(EdgelToolsTest.class.getName());
	
	@Test
	public void testIntervalViewChange(){
		
		Edgel e = new Edgel(
				new double[]{7f, 9f, 11f},
				new double[]{0f, 1f, 0f },
				1f
			);
		
		IntervalIterator srcI = new IntervalIterator(
				new long[]{0,0,0},new long[]{20,20,20} );
		
		IntervalIterator tgtI = new IntervalIterator(
				new long[]{1,1,1},new long[]{20,20,20} );
				
		Edgel eTgtView = EdgelTools.edgelFromInterval(e, srcI, tgtI);
		
		double[] pos = new double[3];
		double[] vpos = new double[3];
		eTgtView.localize( vpos );
		e.localize( pos );
		
		for( int i=0; i<3; i++){
			assertEquals( "pos", pos[i]-1, vpos[i], tol);
			
		}
	}
	
	@Test
	public void testXfm() throws InterruptedException{
		
		Edgel e = new Edgel(
				new double[]{7f, 9f, 11f},
				new double[]{0f, 1f, 0f },
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
		
		Img<FloatType> img = ImgOps.createGradientImgY( new int[]{15, 15, 15},  new FloatType());

		RandomAccess<FloatType> imgRa = img.randomAccess();
		
		int[] patchSize = new int[]{5,5,3};
		Edgel e = new Edgel(
				new double[]{7f, 9f, 11f},
				new double[]{0f, 1f, 0f },
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

	@Test
	public void testZeroXingClose()
	{
		ArrayImgFactory<FloatType> fact = new ArrayImgFactory<FloatType>();
		Img<FloatType> img = fact.create( new int[]{19} , new FloatType());

		float[] dat = new float[]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				2.1e-8f, 0.019f, 0.215f, 4.11f, 17.434f, 32.553f, -17.65f };


		RandomAccess<FloatType> ra = img.randomAccess();
		for( int i=0; i<19; i++){
			ra.setPosition(i, 0);
			ra.get().set(dat[i]);
		}

		double xing = EdgelTools.zeroXing1d(img);

		logger.info("xing: " + xing);
		
	}
	
	@Test
	public void testLaplacianEdge1dClose()
	{
		Img<FloatType> img = ImgOps.createEdgeImg( new int[]{21},
				new double[]{1}, new FloatType(), 1);
		Img<FloatType> lapl = img.factory().create( new int[]{21}, new FloatType());
		Img<FloatType> edge = img.factory().create( new int[]{21}, new FloatType());
		
		img = img.factory().create( new int[]{21}, new FloatType());
		
		img.randomAccess().setPosition(0, 0);
		img.randomAccess().get().setOne();
		System.out.println("\nimg:");
		print(img, 0);
		
		EdgelTools.laplacian( Views.extendBorder( img ), lapl);
		System.out.println("\nlapl:");
		print(lapl, 0);
		
		double edgeX = EdgelTools.zeroXing1d(lapl);
		
		System.out.println("\nedgeX: " + edgeX);
		
		assertEquals( edgeX, 0.5, tol);
	}
	
	public static <T> void print(Img<T> ra, int dim){
		RandomAccess<T> cevRa = ra.randomAccess();
		
		for( int i=0; i < ra.dimension(dim); i++){
			cevRa.setPosition( i , dim );
			System.out.println(" " + cevRa.get());	
		}
	}
	
	@Test
	public void testLaplacianEdge1d()
	{
		Img<FloatType> img = ImgOps.createEdgeImg( new int[]{21},
				new double[]{1}, new FloatType(), 1);
		
		Img<FloatType> lapl = img.factory().create( new int[]{21}, new FloatType());
		Img<FloatType> edge = img.factory().create( new int[]{21}, new FloatType());
		
		
		EdgelTools.laplacian( Views.extendBorder( img ), lapl);
		EdgelTools.localAbsoluteMinNN( 
				Views.extendMirrorDouble( lapl ), 
				Views.interval( edge, Intervals.expand(edge, -1) ), 
				0.0);
		
		RandomAccess<FloatType> ira = img.randomAccess();
		RandomAccess<FloatType> lra = lapl.randomAccess();
		RandomAccess<FloatType> era = edge.randomAccess();
		
		for( int i=0; i< img.dimension(0); i++)
		{
			ira.setPosition(i, 0);
			lra.setPosition(i, 0);
			era.setPosition(i, 0);
			
			
			System.out.println(
					"  " + ira.get() +
					"\t" + lra.get() +
					"\t" + era.get()
			);
			
		}
		
		double zeroXing = EdgelTools.zeroXing1d(lapl);
		System.out.println( "zero x-ing: " + zeroXing);
	}
	
	@Test
	public void testLaplacianEdge3d()
	{
		
		Img<FloatType> img = ImgOps.createEdgeImg( new int[]{15, 15, 15},
				new double[]{1,1,0.1}, new FloatType(), 1);
		
		ImgOps.writeFloat(img, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/edgeImg.tif");
	
		
		ArrayImgFactory<FloatType> factory = new ArrayImgFactory<FloatType>();
		Img<FloatType> lapl = factory.create( new int[]{15, 15, 15}, new FloatType());
		Img<FloatType> edge = factory.create( new int[]{15, 15, 15}, new FloatType());
		
		EdgelTools.laplacian( Views.extendBorder( img ), lapl);
		
		ImgOps.writeFloat(lapl, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/edgeImgLap.tif");
		
//		long[] min = new long[edge.numDimensions()];
//		long[] max = new long[edge.numDimensions()];
//		for ( int d=0; d<edge.numDimensions() - 1; d++ ){
//			min[d] = 1;
//			max[d] = edge.dimension(d) - 1;
//		}
//		min[edge.numDimensions() - 1] = 0;
//		max[edge.numDimensions() - 1] = edge.dimension( edge.numDimensions() - 1);
//		
//		System.out.println("min: " + ArrayUtil.printArray(min));
//		System.out.println("max: " + ArrayUtil.printArray(max));

		
		EdgelTools.localAbsoluteMinNN( 
				Views.extendMirrorDouble( lapl ), 
				Views.interval( edge, Intervals.expand(edge, -1) ), 
				0.000001);
		ImgOps.writeFloat(edge, "/groups/jain/home/bogovicj/projects/crackPatching/toyData/edgeImgLapEdge.tif");
		
	}

}
