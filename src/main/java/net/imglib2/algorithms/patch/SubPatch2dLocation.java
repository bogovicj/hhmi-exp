package net.imglib2.algorithms.patch;

public class SubPatch2dLocation implements Comparable<SubPatch2dLocation> {

	public final int dim;
	public final int xyz;
	
	public final double val;
	
	public SubPatch2dLocation( int dim, int xyz, double val ){
		this.dim = dim;
		this.xyz = xyz;
		this.val = val;
	}

	
	/**
	 * Returns zero if and only if 
	 * this.val == o.val &&
	 * this.dim == o.dim && 
	 * this.xyz == o.xyz
	 * 
	 *    
	 */
	@Override
	public int compareTo(SubPatch2dLocation o) {
		
		if( val > o.val ){
			return  1;
		}else if(val < o.val ){
			return -1;
		}
		
		if( dim > o.dim ){
			return  1;
		}else if( dim < o.dim ){
			return -1;
		}
		
		if( xyz > o.xyz ){
			return  1;
		}else if( xyz < o.xyz ){
			return -1;
		}
		
		return 0;
	}
	
	public SubPatch2dLocation copy(){
		return new SubPatch2dLocation( dim, xyz, val );
	}
	
	
}
