package net.imglib2.realtransform.spline;

import org.junit.Test;
import static org.junit.Assert.*;

import java.util.ArrayList;


import net.imglib2.Cursor;
import net.imglib2.RealLocalizable;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import edu.jhu.ece.iacl.utility.ArrayUtil;

public class ThinPlateSplineKernelTransformTest {
	
	public static double tol = 0.0001;

	public static <T> ArrayList<RealLocalizable> genPtList(Img<T> img){

		ArrayList<RealLocalizable> srcPtList = new ArrayList<RealLocalizable>();
		ArrayList<RealLocalizable> tgtPtList = new ArrayList<RealLocalizable>();

		Cursor<T> front = Views.flatIterable(
				Views.interval( img, Intervals.translate( img, 2, 0 ) ) ).cursor();

		Cursor<T> c = img.localizingCursor();

		while (c.hasNext()) {
			c.next();
			front.next();
			srcPtList.add(c.copyCursor());
			tgtPtList.add(front.copyCursor());
		}

		return srcPtList;
	}
	
//	@Test
//	public void testTranslation(){
//		assertEquals("booyakasha", 50, 50, 0.1);
//	}

	public static ArrayList<RealLocalizable> genPtListStretch(boolean isTgt){
		ArrayList<RealLocalizable> list = new ArrayList<RealLocalizable>();

		double[][] src = new double[][]
				{
				{0,0},
				{0,1},
				{0,2},
				{1,0},
				{1,1},
				{1,2},
				{2,0},
				{2,1},
				{2,2},
				};

		double[][] tgt= new double[][]
				{
				{-0.5, -0.5},
				{-0.5,  1.5},
				{-0.5,  2.0},
				{ 1.5, -0.5},
				{ 1.5,  1.5},
				{ 1.5,  2.0},
				{ 2.0, -0.5},
				{ 2.0,  1.5},
				{ 2.0,  2.0},
				};
		int N = src.length;

		for(int i=0; i<N; i++) {
			if(isTgt){
				list.add( new RL( tgt[i] ));
			}else{
				list.add( new RL( src[i] ));
			}
		}


		return list;
	}

	public static ArrayList<RealLocalizable> genPtListScale(boolean isTgt){
		ArrayList<RealLocalizable> list = new ArrayList<RealLocalizable>();

		for(double x=0; x<3; x++) for(double y=0; y<3; y++){
			if(isTgt){
				list.add(new RL(new double[]{ 2*x , 0.5*y }));
			}else{
				list.add(new RL(new double[]{x,y}));
			}
		}

		return list;
	}

	@Test
	public void testScale() {
		System.out.println("starting");

		ArrayImgFactory<FloatType> factory = new ArrayImgFactory<FloatType>();
		FloatType t = new FloatType();

		Img<FloatType> img = factory.create(new long[] { 3, 3 }, t);

		ArrayList<RealLocalizable> srcPtList = genPtListScale(false);
		ArrayList<RealLocalizable> tgtPtList = genPtListScale(true);

		ThinPlateSplineKernelTransform tps = new ThinPlateSplineKernelTransform( img, srcPtList, tgtPtList);
		tps.computeW();

//		tps.printLandmarks();

		double[] srcPt = new double[]{0.0f,0.0f};
		double[] ptXfm = tps.transformPoint(srcPt);
//		System.out.println("ptXfm: " + ArrayUtil.printArray(srcPt)+ " -> " + ArrayUtil.printArray(ptXfm));
		assertEquals("warp x1", ptXfm[0], 0, tol);
		assertEquals("warp y1", ptXfm[1], 0, tol);

		srcPt = new double[]{0.5f,0.5f};
		ptXfm = tps.transformPoint(srcPt);
//		System.out.println("ptXfm: " + ArrayUtil.printArray(srcPt)+ " -> " + ArrayUtil.printArray(ptXfm));
		assertEquals("warp x2", ptXfm[0], 1.00, tol);
		assertEquals("warp y2", ptXfm[1], 0.25, tol);
		
		srcPt = new double[]{1.0f,1.0f};
		ptXfm = tps.transformPoint(srcPt);
//		System.out.println("ptXfm: " + ArrayUtil.printArray(srcPt)+ " -> " + ArrayUtil.printArray(ptXfm));
		assertEquals("warp x3", ptXfm[0], 2.0, tol);
		assertEquals("warp y3", ptXfm[1], 0.5, tol);
		
	}
	
	@Test
	public void testSimpleWarp() {
		System.out.println("starting");

		ArrayImgFactory<FloatType> factory = new ArrayImgFactory<FloatType>();
		FloatType t = new FloatType();

		Img<FloatType> img = factory.create(new long[] { 3, 3 }, t);

		ArrayList<RealLocalizable> srcPtList = genPtListStretch(false);
		ArrayList<RealLocalizable> tgtPtList = genPtListStretch(true);

		ThinPlateSplineKernelTransform tps = new ThinPlateSplineKernelTransform( img, srcPtList, tgtPtList);
		tps.computeW();

//		tps.printLandmarks();

		double[] srcPt = new double[]{0.0f,0.0f};
		double[] ptXfm = tps.transformPoint(srcPt);
//		System.out.println("ptXfm: " + ArrayUtil.printArray(srcPt)+ " -> " + ArrayUtil.printArray(ptXfm));
		assertEquals("warp x1", ptXfm[0], -0.25, tol);
		assertEquals("warp y1", ptXfm[1], -0.25, tol);

		srcPt = new double[]{0.5f,0.5f};
		ptXfm = tps.transformPoint(srcPt);
//		System.out.println("ptXfm: " + ArrayUtil.printArray(srcPt)+ " -> " + ArrayUtil.printArray(ptXfm));
		assertEquals("warp x2", ptXfm[0], 0.375, tol);
		assertEquals("warp y2", ptXfm[1], 0.375, tol);
		
		srcPt = new double[]{1.0f,1.0f};
		ptXfm = tps.transformPoint(srcPt);
//		System.out.println("ptXfm: " + ArrayUtil.printArray(srcPt)+ " -> " + ArrayUtil.printArray(ptXfm));
		assertEquals("warp x3", ptXfm[0], 1.0, tol);
		assertEquals("warp y3", ptXfm[1], 1.0, tol);
		
		System.out.println("finished");
	}

	public static class RL implements RealLocalizable{
		double[] d;
		public RL(double[] d){
			this.d = d;
		}

		public int numDimensions() {
			return d.length;
		}

		public double getDoublePosition(int i) {
			return d[i];
		}

		public float getFloatPosition(int i) {
			return (float)d[i];
		}

		public void localize(float[] loc) {
			for(int i=0; i<loc.length; i++){
				loc[i] = getFloatPosition(i);
			}
		}

		public void localize(double[] loc) {
			for(int i=0; i<loc.length; i++){
				loc[i] = getDoublePosition(i);
			}
		}
	}
}
