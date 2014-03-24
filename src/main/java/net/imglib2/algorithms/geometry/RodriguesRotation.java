package net.imglib2.algorithms.geometry;

import org.ejml.simple.SimpleMatrix;

import edu.jhu.ece.iacl.utility.ArrayUtil;
import net.imglib2.realtransform.AffineTransform3D;

public class RodriguesRotation {
	
	
	/**
	 * 
	 * @param src source vector
	 * @param dest destination vector
	 * @return a transform that maps the src vector to the dest vector
	 */
	public static AffineTransform3D rotation(double[] src, double[] dest){
		
		if( src.length != 3 || dest.length !=3 ){
			return null;
		}
		
		double[] srcNorm = ArrayUtil.normalizeLength(src);
		double[] destNorm = ArrayUtil.normalizeLength(dest);
		
		double[] cross = crossProduct3d( srcNorm, destNorm );
		
		SimpleMatrix srcM = new SimpleMatrix( new double[][]{ srcNorm  } );
		SimpleMatrix dstM = new SimpleMatrix( new double[][]{ destNorm } );
		SimpleMatrix crsM = new SimpleMatrix( new double[][]{ cross    } );	
		SimpleMatrix I = SimpleMatrix.identity(3);
		
		double costheta = srcM.mult(dstM.transpose()).get(0);
		
		SimpleMatrix k = new SimpleMatrix(3,3);
		k.set(0, 1, -crsM.get( 0, 2 ));
		k.set(0, 2,  crsM.get( 0, 1 ));
		
		k.set(1, 0,  crsM.get( 0, 2 ));
		k.set(1, 2, -crsM.get( 0, 1 ));
		
		k.set(2, 0, -crsM.get( 0, 1 ));
		k.set(2, 1,  crsM.get( 0, 0 ));
		
		SimpleMatrix R = 
				I.scale(costheta).
				plus(k).
				plus( k.mult(k.transpose()).scale( 1 - costheta ).scale(1/sumSquared(k)) ); 
		
		System.out.println("R: " + R);
		
		AffineTransform3D out = new AffineTransform3D();
		out.set( R.get(0, 0), R.get(0, 1), R.get(0, 2), 0, 
				 R.get(1, 0), R.get(1, 1), R.get(1, 2), 0,
				 R.get(2, 0), R.get(2, 1), R.get(2, 2), 0 );
		
		return out;
	}
	
	public static double sumSquared(SimpleMatrix mtx){
		
		double out = 0;
		
		int N = mtx.numRows()*mtx.numCols();
		for(int i=0; i<N; i++){
			out += mtx.get(i)* mtx.get(i);
		}
		
		return out;
	}
	
	public static double[] crossProduct3d(double[] x, double[] y){
		
		if( x.length != 3 || y.length !=3 ){
			return null;
		}
		double[] z = new double[3];
		
		z[0] = x[1]*y[2] - x[2]*y[1];
		z[1] = x[2]*y[0] - x[0]*y[2];
		z[2] = x[0]*y[1] - x[1]*y[0];
		
		return z;
	}

	public static void main(String[] args){
		
		double[] u = new double[]{ 0.245, -0.563, -0.055 };
		double[] v = new double[]{ 0.3,   -0.5,    0.2   };
		
		AffineTransform3D mtx = rotation(u,v);
		
		System.out.println("mtx: " + mtx);
		
		double[] uXfm = new double[3];
		mtx.apply(u, uXfm);
		
		System.out.println("u   : " + ArrayUtil.printArray(u));
		System.out.println("uXfm: " + ArrayUtil.printArray(uXfm));
		
	}
}
