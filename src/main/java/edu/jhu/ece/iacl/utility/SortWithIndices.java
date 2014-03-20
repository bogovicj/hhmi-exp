package edu.jhu.ece.iacl.utility;

import java.util.Arrays;


/**
 * 
 * Taken from stackoverflow 
 * http://stackoverflow.com/questions/951848/java-array-sort-quick-way-to-get-a-sorted-list-of-indices-of-an-array
 * 
 * Usage:
 * double[] x = ... // the array to sort
 * double[] i = SortWithIndices.getIndexArray(x.length); // generate an array with [0, 1, ..., x.length]
 * SortWithIndices.quicksort(x,i);  // sorts x in-place
 *
 *  @author kd304 (stackoverflow)
 *	@author John Bogovic
 */
public class SortWithIndices {
	
	
	public static void quicksort(float[] main, int[] origindex) {
	    quicksort(main, origindex, 0, origindex.length - 1);
	}
	public static int[] getBIndexArray(int N){
		int[] idxs = new int[N];
		for(int i=0; i<N; i++){
			idxs[i]=N-1-i;
		}
		return idxs;
	}
	
	public static int[] getIndexArray(int N){
		int[] idxs = new int[N];
		for(int i=0; i<N; i++){
			idxs[i]=i;
		}
		return idxs;
	}
	
	public static int[] reorder(int[] in, int[] idx){
		int[] out = new int[idx.length];
		for(int i=0; i<idx.length; i++){
			out[i]=in[idx[i]];
		}
		return out;
	}
	public static float[] reorder(float[] in, int[] idx){
		float[] out = new float[idx.length];
		for(int i=0; i<idx.length; i++){
			out[i]=in[idx[i]];
		}
		return out;
	}
	public static double[] reorder(double[] in, int[] idx){
		double[] out = new double[idx.length];
		for(int i=0; i<idx.length; i++){
			out[i]=in[idx[i]];
		}
		return out;
	}
	
	public static int[] reorderReverse(int[] in, int[] idx){
		int[] out = new int[idx.length];
		int j = 0;
		for(int i=idx.length-1; i>=0; i--){
			out[j]=in[idx[i]];
			j++;
		}
		return out;
	}
	public static float[] reorderReverse(float[] in, int[] idx){
		float[] out = new float[idx.length];
		int j = 0;
		for(int i=idx.length-1; i>=0; i--){
			out[j]=in[idx[i]];
			j++;
		}
		return out;
	}
	public static double[] reorderReverse(double[] in, int[] idx){
		double[] out = new double[idx.length];
		int j = 0;
		for(int i=idx.length-1; i>=0; i--){
			out[j]=in[idx[i]];
			j++;
		}
		return out;
	}
	
	// quicksort a[left] to a[right]
	public static void quicksort(float[] a, int[] origindex, int left, int right) {
	    if (right <= left) return;
	    int i = partition(a, origindex, left, right);
	    quicksort(a, origindex, left, i-1);
	    quicksort(a, origindex, i+1, right);
	}

	// partition a[left] to a[right], assumes left < right
	private static int partition(float[] a, int[] origindex, 
	int left, int right) {
	    int i = left - 1;
	    int j = right;
	    while (true) {
	        while (less(a[++i], a[right]))      // find item on left to swap
	            ;                               // a[right] acts as sentinel
	        while (less(a[right], a[--j]))      // find item on right to swap
	            if (j == left) break;           // don't go out-of-bounds
	        if (i >= j) break;                  // check if pointers cross
	        exch(a, origindex, i, j);               // swap two elements into place
	    }
	    exch(a, origindex, i, right);               // swap with partition element
	    return i;
	}

	// is x < y ?
	private static boolean less(float x, float y) {
	    return (x < y);
	}

	// exchange a[i] and a[j]
	private static void exch(float[] a, int[] index, int i, int j) {
		float swap = a[i];
	    a[i] = a[j];
	    a[j] = swap;
	    int b = index[i];
	    index[i] = index[j];
	    index[j] = b;
	}
	
	
	
	
	public static void quicksort(int[] main, int[] origindex) {
	    quicksort(main, origindex, 0, origindex.length - 1);
	}

	// quicksort a[left] to a[right]
	public static void quicksort(int[] a, int[] origindex, int left, int right) {
	    if (right <= left) return;
	    int i = partition(a, origindex, left, right);
	    quicksort(a, origindex, left, i-1);
	    quicksort(a, origindex, i+1, right);
	}

	// partition a[left] to a[right], assumes left < right
	private static int partition(int[] a, int[] origindex, 
	int left, int right) {
	    int i = left - 1;
	    int j = right;
	    while (true) {
	        while (less(a[++i], a[right]))      // find item on left to swap
	            ;                               // a[right] acts as sentinel
	        while (less(a[right], a[--j]))      // find item on right to swap
	            if (j == left) break;           // don't go out-of-bounds
	        if (i >= j) break;                  // check if pointers cross
	        exch(a, origindex, i, j);               // swap two elements into place
	    }
	    exch(a, origindex, i, right);               // swap with partition element
	    return i;
	}

	// is x < y ?
	private static boolean less(int x, int y) {
	    return (x < y);
	}

	
	// exchange a[i] and a[j]
	private static void exch(int[] a,int[] origindex, int i, int j) {
		int swap = a[i];
	    a[i] = a[j];
	    a[j] = swap;
	    
	    int b = origindex[i];
	    origindex[i] = origindex[j];
	    origindex[j] = b;
	    
	} 
	
	
	public static void quicksort(double[] main, int[] origindex) {
	    quicksort(main, origindex, 0, origindex.length - 1);
	}

	// quicksort a[left] to a[right]
	public static void quicksort(double[] a, int[] origindex, int left, int right) {
	    if (right <= left) return;
	    int i = partition(a, origindex, left, right);
	    quicksort(a, origindex, left, i-1);
	    quicksort(a, origindex, i+1, right);
	}

	// partition a[left] to a[right], assumes left < right
	private static int partition(double[] a, int[] origindex, 
	int left, int right) {
	    int i = left - 1;
	    int j = right;
	    while (true) {
	        while (less(a[++i], a[right]))      // find item on left to swap
	            ;                               // a[right] acts as sentinel
	        while (less(a[right], a[--j]))      // find item on right to swap
	            if (j == left) break;           // don't go out-of-bounds
	        if (i >= j) break;                  // check if pointers cross
	        exch(a, origindex, i, j);               // swap two elements into place
	    }
	    exch(a, origindex, i, right);               // swap with partition element
	    return i;
	}

	// is x < y ?
	private static boolean less(double x, double y) {
	    return (x < y);
	}
	// is x < y ?
	private static boolean greater(double x, double y) {
		return (x > y);
	}
	
	// exchange a[i] and a[j]
	private static void exch(double[] a,int[] origindex, int i, int j) {
		double swap = a[i];
	    a[i] = a[j];
	    a[j] = swap;
	    
	    int b = origindex[i];
	    origindex[i] = origindex[j];
	    origindex[j] = b;
	    
	} 
	
	public static void main(String[] args){
		
		int[] x = new int[]{ 5, 9, 1, 21, 104, 2};
		int[] idxs = getIndexArray(x.length);
		
		System.out.println("x: " + Arrays.toString(x));
		System.out.println("idxs: " + Arrays.toString(idxs));
		
		quicksort(x,idxs);
		System.out.println("x: " + Arrays.toString(x));
		System.out.println("idxs: " + Arrays.toString(idxs));
		
		
	}


}
