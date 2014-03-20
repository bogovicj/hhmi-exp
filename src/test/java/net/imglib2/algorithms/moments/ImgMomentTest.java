package net.imglib2.algorithms.moments;

import static org.junit.Assert.*;

import org.ejml.data.DenseMatrix64F;
import org.junit.Test;

import java.util.Arrays;

import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.type.numeric.real.FloatType;
import edu.jhu.ece.iacl.utility.ArrayUtil;

public class ImgMomentTest {

//	@Test 
//	public void testBlah(){
//
//		int[] permutation = new int[]{0,1};
//
//		int N = 2;
//
//		DenseMatrix64F out = new DenseMatrix64F(2,2);
//		out.set(0, 0, 0);
//		out.set(0, 1, 1);
//		out.set(1, 0, 2);
//		out.set(1, 1, 3);
//		
//		for(int i=0; i<out.numRows; i++){
//			for(int j=0; j<out.numCols; j++){
//				
//				DenseMatrix64F evec = in.getEigenVector(permutation[out.numCols - j - 1]);
//				out.set(i,j,in.getEigenVector(permutation[out.numCols - j - 1]).get(i));
//				
//			}
//		}
//	}
	
	@Test
	public void testOrientationsSimple2d(){

		FloatType t = new FloatType(); 
		ArrayImg<FloatType,FloatType> img = new ArrayImg<FloatType,FloatType>(t, new long[]{5,5}, 1); 	

		Img<FloatType> a = img.factory().create(new long[]{5,5}, t);
		Img<FloatType> b = img.factory().create(new long[]{5,5}, t);
		Img<FloatType> c = img.factory().create(new long[]{5,5}, t);
		//System.out.println("a " + a);

		RandomAccess<FloatType> ra = a.randomAccess();
		RandomAccess<FloatType> rb = b.randomAccess();
		RandomAccess<FloatType> rc = c.randomAccess();
				
	 	for (int i=0; i<5; i++){
			ra.setPosition(new long[]{i,2});
			rb.setPosition(new long[]{2,i});

			ra.get().setReal(1);
			rb.get().setReal(1);
		}

		ImgMoment mom_a = new ImgMoment();
		ImgMoment mom_b = new ImgMoment();

		double[] or_a = mom_a.orientation(a);
		double[] or_b = mom_b.orientation(b);

		System.out.println("or_a: " + ArrayUtil.printArray(or_a));
		System.out.println("or_b: " + ArrayUtil.printArray(or_b));

		assertTrue("x cov 2d", Arrays.equals(or_a, new double[]{2,0,0,0}));	
		assertTrue("y cov2d", Arrays.equals(or_b, new double[]{0,0,0,2}));	

		mom_a.orientationEvalsEvecs(or_a);
		mom_b.orientationEvalsEvecs(or_b);

//		System.out.println("eval_a: \n" + ArrayUtil.printArray(mom_a.getEvals()));
//		System.out.println("eval_b: \n" + ArrayUtil.printArray(mom_b.getEvals()));
//
		System.out.println("evec_a: \n" + ArrayUtil.printArray(ArrayUtil.reshape2D(mom_a.getEvecs(),2,2,false)));
		System.out.println("evec_b: \n" + ArrayUtil.printArray(ArrayUtil.reshape2D(mom_b.getEvecs(),2,2,false)));

	}
	
	@Test
	public void testOrientationsSimple3d(){

		FloatType t = new FloatType(); 
		ArrayImg<FloatType,FloatType> img = new ArrayImg<FloatType,FloatType>(t, new long[]{5,5,5}, 1); 	

		Img<FloatType> a = img.factory().create(new long[]{5,5,5}, t);
		Img<FloatType> b = img.factory().create(new long[]{5,5,5}, t);
		Img<FloatType> c = img.factory().create(new long[]{5,5,5}, t);
		//System.out.println("a " + a);

		RandomAccess<FloatType> ra = a.randomAccess();
		RandomAccess<FloatType> rb = b.randomAccess();
		RandomAccess<FloatType> rc = c.randomAccess();
				
	 	for (int i=0; i<5; i++){
			ra.setPosition(new long[]{i,2,2});
			rb.setPosition(new long[]{2,i,2});
			rc.setPosition(new long[]{2,2,i});

			ra.get().setReal(1);
			rb.get().setReal(1);
			rc.get().setReal(1);
		}

		ImgMoment mom_a = new ImgMoment();
		ImgMoment mom_b = new ImgMoment();
		ImgMoment mom_c = new ImgMoment();

		double[] or_a = mom_a.orientation(a);
		double[] or_b = mom_b.orientation(b);
		double[] or_c = mom_c.orientation(c);


		//System.out.println("or_a: " + ArrayUtil.printArray(or_a));
		//System.out.println("or_b: " + ArrayUtil.printArray(or_b));
		//System.out.println("or_c: " + ArrayUtil.printArray(or_c));

		assertTrue("x cov 3d", Arrays.equals(or_a, new double[]{2,0,0,0,0,0,0,0,0}));	
		assertTrue("y cov 3d", Arrays.equals(or_b, new double[]{0,0,0,0,2,0,0,0,0}));	
		assertTrue("z cov 3d", Arrays.equals(or_c, new double[]{0,0,0,0,0,0,0,0,2}));	

		mom_a.orientationEvalsEvecs(or_a);
		mom_b.orientationEvalsEvecs(or_b);
		mom_c.orientationEvalsEvecs(or_c);

		System.out.println("eval_a: " + ArrayUtil.printArray(mom_a.getEvals()));
		System.out.println("eval_b: " + ArrayUtil.printArray(mom_b.getEvals()));
		System.out.println("eval_c: " + ArrayUtil.printArray(mom_c.getEvals()));

		System.out.println("evec_a:\n" + ArrayUtil.printArray(ArrayUtil.reshape2D(mom_a.getEvecs(),3,3,false)));
		System.out.println("evec_b:\n" + ArrayUtil.printArray(ArrayUtil.reshape2D(mom_b.getEvecs(),3,3,false)));
		System.out.println("evec_c:\n" + ArrayUtil.printArray(ArrayUtil.reshape2D(mom_c.getEvecs(),3,3,false)));

	}

}
