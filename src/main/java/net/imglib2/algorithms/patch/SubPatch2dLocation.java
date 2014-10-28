package net.imglib2.algorithms.patch;

import java.io.Serializable;

public class SubPatch2dLocation implements Comparable<SubPatch2dLocation>, Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 271503445605322074L;
	public final int dim;
	public final int xyz;
	
	public final int idx;
	public final double val;
	
	public SubPatch2dLocation( int dim, int xyz, int idx, double val ){
		this.dim = dim;
		this.xyz = xyz;
		this.idx = idx;
		this.val = val;
	}

	
	/**
	 * Returns zero if and only if 
	 * this.val == o.val &&
	 * this.dim == o.dim &&
	 * this.idx == o.idx && 
	 * this.xyz == o.xyz
	 * 
	 *    
	 */
	@Override
	public int compareTo(SubPatch2dLocation o) {
		
		// sort by 
		if( val > o.val ){
			return  1;
		}else if(val < o.val ){
			return -1;
		}

		
		if ( dim != o.dim ||
			 xyz != o.xyz ||
			 idx != o.idx )
		{
			if( Math.random() > 0.5 ){
				return 1;
			}else{
				return -1;
			}
		}
		
		
		/*
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
		*/
		
		return 0;
	}
	
	public SubPatch2dLocation copy(){
		return new SubPatch2dLocation( dim, xyz, idx, val );
	}
	
	public String toString(){
		return String.format("SPL dim(%d) xyz(%d) idx(%d) val(%f)", this.dim, this.xyz, this.idx, this.val );
	}
	
	
}
