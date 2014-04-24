package net.imglib2.algorithms.patch;

import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import edu.jhu.ece.iacl.utility.ArrayUtil;

public class PatchTools {


	/**
	 * Returns the integer coordinate of the midpoint of a patch
	 * with the specified size.  Is only accurate for patches
	 * with an odd size in every dimension.
	 * 
	 * @param patchSize size of the patch
	 * @return the midpoint coordinate
	 */
	public static int[] patchSizeToMidpt(int[] patchSize){ // determine translation

		int[] midPt = ArrayUtil.clone(patchSize);
		ArrayUtil.addInPlace(midPt, -1);
		ArrayUtil.divide(midPt, 2);

		return midPt;
	}
	
	public static <T> void copyViewTo( RandomAccessible<T> src, RandomAccessibleInterval<T> dest)
	{
		
	}
	
		
}
