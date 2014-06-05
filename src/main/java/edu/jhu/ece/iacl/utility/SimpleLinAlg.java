package edu.jhu.ece.iacl.utility;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

public class SimpleLinAlg {

	Logger logger = LogManager.getLogger(SimpleLinAlg.class.getName());
	
	/**
	 * returns a 2d vector orthogonal to the input.
	 * @param vec
	 * @return
	 */
	public static float[] orth2d( float[] vec ){
		if( vec.length != 2){
			return null;
		}
		float axSqr = vec[0]*vec[0];
		float aySqr = vec[1]*vec[1];
		
		float[] out = new float[2];
		out[1] = (float) Math.sqrt ( axSqr / ( axSqr + aySqr) );
		out[0] = - vec[1] * out[1] / vec[0];
 		
		return out;
	}
	
	/**
	 * returns a 2d vector orthogonal to the input.
	 * @param vec
	 * @return
	 */
	public static double[] orth2d( double[] vec ){
		if( vec.length != 2){
			return null;
		}
		double axSqr = vec[0]*vec[0];
		double aySqr = vec[1]*vec[1];
		
		double[] out = new double[2];
		out[1] = Math.sqrt ( axSqr / ( axSqr + aySqr) );
		out[0] = - vec[1] * out[1] / vec[0];
 			
		return out;
	}
	
	public static void main( String[] args ){
		float[] u = new float[]{2f, 1f};
		float[] v = orth2d( u );
		System.out.println(" v: " + ArrayUtil.printArray( v ));
	}
}
	