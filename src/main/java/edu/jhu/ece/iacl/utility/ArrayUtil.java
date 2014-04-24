package edu.jhu.ece.iacl.utility;

import java.util.List;
import java.util.Arrays;


/**
*
*  This class computes various basic functionality for dealing with arrays.
*	<p> 
*	Includes indexing, some arithmetic, normalization, reshaping, conversion to Jama Matrix or JIST ImageData
*
*	@version    July 20, 2011
*	@author     John Bogovic
*/

public class ArrayUtil {
//	private static int testx = -1;
//	private static int testy = -1;
//	private static int testz = -1;

	private double[] minbounds;
	private double[] maxbounds;
	private double[] step;
	private double[] currentVector;
	private int movingComp;
	private int numComps;
	
	public ArrayUtil(){}
	
	public void initVectorSweep(double[] minbounds, double[] maxbounds, double[] step){
		this.minbounds = minbounds;
		this.maxbounds = maxbounds;
		this.step = step;
		if(minbounds.length!=maxbounds.length){
			System.out.println("Warning - vectors must be of the same length");
		}
		if(minbounds.length!=step.length){
			System.out.println("Warning - vectors must be of the same length");
		}
		numComps = minbounds.length;
		System.out.println("numcomps: " + numComps);
	}
	
	public static void testDouble(double a){
		System.out.println(a);
	}
	
	public static void testInt(int a){
		System.out.println(a);
	}
	
	public static void test2dFloatArray(float[][] a){
		System.out.println(a[0][0]);
		System.out.println(a[0][0]);
	}
	
	public double[] getNextVector(){
		// init current Vector if needed
		if(currentVector==null){
			currentVector=minbounds.clone();
			movingComp = numComps-1;
			return currentVector;
		}
		currentVector[movingComp]=currentVector[movingComp]+step[movingComp];
		
		while(currentVector[movingComp]>maxbounds[movingComp]){
			currentVector[movingComp]=minbounds[movingComp];
			movingComp--;
			
			if(movingComp<0){
				return null;
			}else{
				currentVector[movingComp]=currentVector[movingComp]+step[movingComp];
			}
		}
		
		movingComp = numComps-1;
		return currentVector;
	}



	public static final int indexOf(int x, int[] list){
		for(int i=0; i<list.length; i++){
			if(x==list[i])
				return i;
		}
		return -1;
	}
	public static final int indexOf(float x, float[] list){
		for(int i=0; i<list.length; i++){
			if(x==list[i])
				return i;
		}
		return -1;
	}
	public static final int indexOf(double x, double[] list){
		for(int i=0; i<list.length; i++){
			if(x==list[i])
				return i;
		}
		return -1;
	}
	public static final int indexOf(String x, String[] list){
		for(int i=0; i<list.length; i++){
			if(x.equals(list[i]))
				return i;
		}
		return -1;
	}
	
	/**
	 *  Returns all the indices of x in list
	 */
	public static final int[] find(int x, int[] list){
		int n = numberOf(x,list);
		if(n==0){ return null; }
		
		int[] idxs = new int[n];
		
		int k = 0;
		for(int i=0; i<list.length; i++){
			if(x==list[i]){
				idxs[k]=i;
				k++;
				
				if(k==n){ // finish early if we can
					return idxs;
				}
			}
		}
		
		return idxs;
	}
	
	public static final int[] find(boolean x, boolean[] list){
		int n = numberOf(x,list);
		if(n==0){ return null; }
		
		int[] idxs = new int[n];
		
		int k = 0;
		for(int i=0; i<list.length; i++){
			if(x==list[i]){
				idxs[k]=i;
				k++;
				
				if(k==n){ // finish early if we can
					return idxs;
				}
			}
		}
		
		return idxs;
	}
	
	public static final int[] find(double x, double[] list){
		int n = numberOf(x,list);
		if(n==0){ return null; }
		
		int[] idxs = new int[n];
		
		int k = 0;
		for(int i=0; i<list.length; i++){
			if(x==list[i]){
				idxs[k]=i;
				k++;
				
				if(k==n){ // finish early if we can
					return idxs;
				}
			}
		}
		
		return idxs;
	}
	
	public static final int[] find(String x, String[] list){
		int n = numberOf(x,list);
		if(n==0){ return null; }
		
		int[] idxs = new int[n];
		
		int k = 0;
		for(int i=0; i<list.length; i++){
			if(x.equals(list[i])){
				idxs[k]=i;
				k++;
				
				if(k==n){ // finish early if we can
					return idxs;
				}
			}
		}
		
		return idxs;
	}
	
	public static final int numberOf(boolean x, boolean[] list){
		int n = 0;
		for(int i=0; i<list.length; i++){
			if(x==list[i])
				n++;
		}
		return n;
	}
	
	public static final int numberOf(int x, int[] list){
		int n = 0;
		for(int i=0; i<list.length; i++){
			if(x==list[i])
				n++;
		}
		return n;
	}
	public static final int numberOf(float x, float[] list){
		int n = 0;
		for(int i=0; i<list.length; i++){
			if(x==list[i])
				n++;
		}
		return n;
	}
	public static final int numberOf(double x, double[] list){
		int n = 0;
		for(int i=0; i<list.length; i++){
			if(x==list[i])
				n++;
		}
		return n;
	}
	public static final int numberOf(String x, String[] list){
		int n = 0;
		for(int i=0; i<list.length; i++){
			if(x.equals(list[i]))
				n++;
		}
		return n;
	}
	
	public static final void replace(byte find, byte replace, byte[] list){
		for(int i=0; i<list.length; i++){
			if(list[i]==find)
				list[i]=replace;
		}
	}
	
	public static final int[] sortWithInds(int[] x, int[] idx){
		int[] out = x.clone();
		Arrays.sort(out);
		//populate the index array
		for(int i=0; i<out.length; i++) 
			idx[i]=indexOf(out[i],x);
		
		return out;
	}
	public static final float[] sortWithInds(float[] x, int[] idx){
		float[] out = x.clone();
		Arrays.sort(out);
		//populate the index array
		for(int i=0; i<out.length; i++) 
			idx[i]=indexOf(out[i],x);
		
		return out;
	}
	
	public static final double[] sortWithInds(double[] x, int[] idx, int[] rev){
		double[] out = x.clone();
		Arrays.sort(out);
		//populate the index array
		for(int i=0; i<out.length; i++) 
			idx[i]=indexOf(out[i],x);
		
		return out;
	}
	
	public static final int[] sortWithInds(int[] x, int[] idx, int[] rev){
		int[] out = x.clone();
		Arrays.sort(out);
		//populate the index array
		for(int i=0; i<out.length; i++) 
			idx[i]=indexOf(out[i],x);
		
		return out;
	}
	public static final float[] sortWithInds(float[] x, int[] idx, int[] rev){
		float[] out = x.clone();
		Arrays.sort(out);
		//populate the index array
		for(int i=0; i<out.length; i++) 
			idx[i]=indexOf(out[i],x);
		
		return out;
	}
	public static final double[] sortWithInds(double[] x, int[] idx){
		double[] out = x.clone();
		Arrays.sort(out);
		//populate the index array
		for(int i=0; i<out.length; i++) 
			idx[i]=indexOf(out[i],x);
		
		return out;
	}
    public static float[] reorder(float[] orig, int[] inds){
        float[] out = new float[orig.length];
        for(int i=0; i<orig.length; i++){
            out[i] = orig[inds[i]];
        }
        return out;
    }
    public static int[] reverseInds(int[] inds){
        int N = inds.length;
        int[] revinds = new int[N];
        for(int i=0; i<N; i++){
            revinds[i] = inds[N-inds[i]-1];
        }
        return revinds;
    }
	
	public static final int[] permute(int[] in, int[] idx){
		int[] out = new int[in.length];
		for(int i=0; i<in.length; i++)
			out[i]=in[idx[i]];
		
		return out;
	}
	public static final float[] permute(float[] in, int[] idx){
		float[] out = new float[in.length];
		for(int i=0; i<in.length; i++)
			out[i]=in[idx[i]];
		
		return out;
	}
	public static final double[] permute(double[] in, int[] idx){
		double[] out = new double[in.length];
		for(int i=0; i<in.length; i++)
			out[i]=in[idx[i]];
		
		return out;
	}
	
	public static final float min(float[] in, int[] idx){
		float min = in[0];
		for(int i=1; i<in.length; i++){
			if(in[i]<min){ min=in[i]; }
		}
		return min;
	}
	public static final float max(float[] in, int[] idx){
		float max = in[0];
		for(int i=1; i<in.length; i++){
			if(in[i]>max){ max=in[i]; }
		}
		return max;
	}

	/**
	 * 
	 * @param x 
	 * @param list
	 * @return True if x is contained in the list
 	 */
	public static final boolean contains(int[] list, int x){
		return (ArrayUtil.indexOf(x,list)>-1);
	}
	/**
	 * 
	 * @param x 
	 * @param list
	 * @return True if x is contained in the list
 	 */
	public static final boolean contains(float[] list, float x){
		return (ArrayUtil.indexOf(x,list)>-1);
	}
	
	/**
	 * 
	 * @param x 
	 * @param list
	 * @return True if x is contained in the list
 	 */
	public static final boolean contains(String[] list, String x){
		return (ArrayUtil.indexOf(x,list)>-1);
	}
	
	/**
	 * If list contains x, returns a new array without x.
	 * Otherwise returns list (not a clone)
	 * @param list
	 * @param el
	 * @return
	 */
	public static final int[] removeElement(int[] list, int x){
		if(ArrayUtil.contains(list,x)){
			int[] out = new int[list.length-1];
			int j=0;
			for(int i=0; i<list.length; i++){
				if(list[i]!=x){
					out[j]=list[i];
					j++;
				}
			}
			return out;
		}else{
			return list;
		}
	}
	
	public static final float[] removeIndex(float[] in, int ind){
		float[] out = new float[in.length-1];
		int j=0;
		for(int i=0; i<in.length; i++){
			if(i!=ind){
				out[j]=in[i];
				j++;
			}
		}
		return out;
	}
	public static final int[] removeIndex(int[] in, int ind){
		int[] out = new int[in.length-1];
		int j=0;
		for(int i=0; i<in.length; i++){
			if(i!=ind){
				out[j]=in[i];
				j++;
			}
		}
		return out;
	}
	public static final void fill(byte[] in, byte val){
		for(int i=0; i<in.length; i++){
			in[i] = val;
		}
	}
	public static final void fill(int[] in, int val){
		for(int i=0; i<in.length; i++){
			in[i] = val;
		}
	}
	public static final void fill(float[] in, float val){
		for(int i=0; i<in.length; i++){
			in[i] = val;
		}
	}
	public static final void fill(float[][] in, float val){
		for(int i=0; i<in.length; i++)for(int j=0; j<in[0].length; j++){
			in[i][j] = val;
		}
	}
	public static final void fill(boolean[] in, boolean val){
		for(int i=0; i<in.length; i++){
			in[i] = val;
		}
	}
	public static final void fill(double[] in, double val){
		for(int i=0; i<in.length; i++){
			in[i] = val;
		}
	}
	public static final void fill(boolean[][] in, boolean val){
		for(int i=0; i<in.length; i++)for(int j=0; j<in[0].length; j++){
			in[i][j] = val;
		}
	}
	public static final void fill(double[][] in, double val){
		for(int i=0; i<in.length; i++)for(int j=0; j<in[0].length; j++){
			in[i][j] = val;
		}
	}
	public static final void fill(int[][] in, int val){
		for(int i=0; i<in.length; i++)for(int j=0; j<in[0].length; j++){
			in[i][j] = val;
		}
	}
	public static final void fill(boolean[][][] in, boolean val){
		for(int i=0; i<in.length; i++)for(int j=0; j<in[0].length; j++)for(int k=0; k<in[0][0].length; k++){
			in[i][j][k] = val;
		}
	}
	public static final void fill(float[][][] in, float val){
		for(int i=0; i<in.length; i++)for(int j=0; j<in[0].length; j++)for(int k=0; k<in[0][0].length; k++){
			in[i][j][k] = val;
		}
	}
	public static final void fill(double[][][] in, double val){
		for(int i=0; i<in.length; i++)for(int j=0; j<in[0].length; j++)for(int k=0; k<in[0][0].length; k++){
			in[i][j][k] = val;
		}
	}
	public static final void fill(int[][][] in, int val){
		for(int i=0; i<in.length; i++)for(int j=0; j<in[0].length; j++)for(int k=0; k<in[0][0].length; k++){
			in[i][j][k] = val;
		}
	}
	public static final void fill(float[][][][] in, float val){
		for(int i=0; i<in.length; i++)for(int j=0; j<in[0].length; j++) for(int k=0; k<in[0][0].length; k++) for(int l=0; l<in[0][0][0].length; l++){
			in[i][j][k][0] = val;
		}
	}
	
	/**
	 * Inserts all the elements of a source array into the destination array
	 * where the starting index of 
	 * 
	 * @param dest destination array
	 * @param src source array
	 * @param i starting row position
	 * @param j starting column position
	 */
	public static void insert( double[][] dest, double[][] src, int i, int j)
	{
		int ii=0, jj=0;
		for( int r=i; r<i+src.length; r++)
		{
			for( int c=i; c<j+src[0].length; c++)
			{
				dest[r][c] = src[ii][jj];
				jj++;
			}
			ii++;
		}
	}
	
	public static float[] concatenate(List<float[]> list){
		int totlen = 0;
		for(int i=0; i<list.size(); i++){
			totlen += list.get(i).length;
		}
		float[] out = new float[totlen];
		int k = 0;
		for(int i=0; i<list.size(); i++){
			float[] a = list.get(i);
			for(int j=0; j<a.length; j++){
				out[k]=a[j];
				k++;
			}
		}
		return out;
	}
	
	public static float[][] concatenate(float[][] a, float[][] b, int dim){
		float[][] out = null;
		if(dim==0){
			if(a[0].length!=b[0].length){
				System.out.println("number of columns of inputs are not the identical");
				return null;
			}
			out = new float[a.length + b.length][a[0].length];
			
			int k = 0;
			for(int i=0; i<a.length; i++){
				for(int j=0; j<a[0].length; j++){
					out[i][j]=a[i][j];
				}
				k++;
			}
			for(int i=0; i<b.length; i++){
				for(int j=0; j<a[0].length; j++){
					out[k][j]=b[i][j];
				}
				k++;
			}
			
		}else{
			if(a.length!=b.length){
				System.out.println("number of rows of inputs are not the identical");
				return null;
			}
			out = new float[a.length][a[0].length+ b[0].length];
			
			int k = 0;
			for(int j=0; j<a[0].length; j++){
				for(int i=0; i<a.length; i++){
					out[i][j]=a[i][j];
				}
				k++;
			}
			for(int j=0; j<b[0].length; j++){
				for(int i=0; i<b.length; i++){
					out[i][k]=b[i][j];
				}
				k++;
			}
			
		}
		return out;
	}
	
	public static double[][] concatenate(List<double[][]> alist, int dim){
		double[][] out = null;
		if(dim==0){ // concat rows
			int N = 0;
			int M = 0;
			for(double[][] l : alist){
				N+=l.length;
				if(l[0].length>M){
					M = l[0].length;
				}
			}
			out = new double[N][M];
			
			int k = 0;
			for(double[][] l : alist){
				for(int i=0; i<l.length; i++)for(int j=0; j<l[i].length; j++){
					out[k][j]=l[i][j];
					k++;
				}
			}
			
		}else{ // concat cols
			int N = 0;
			int M = 0;
			for(double[][] l : alist){
				N+=l[0].length;
				if(l.length>M){
					M = l.length;
				}
			}
			out = new double[M][N];
			int k = 0;
			for(double[][] l : alist){
				for(int i=0; i<l.length; i++)for(int j=0; j<l[i].length; j++){
					out[i][k]=l[i][j];
					k++;
				}
			}
		}
		
		return out;
	}
	
	public static final int[] clone(int[] in){
		int nx = in.length;
		int[]out = new int[nx];
		for(int i=0; i<nx; i++){
			out[i] = in[i];
		}
		return out;
	}
	public static final float[] clone(float[] in){
		int nx = in.length;
		float[]out = new float[nx];
		for(int i=0; i<nx; i++){
			out[i] = in[i];
		}
		return out;
	}
	public static final byte[] clone(byte[] in){
		int nx = in.length;
		byte[]out = new byte[nx];
		for(int i=0; i<nx; i++){
			out[i] = in[i];
		}
		return out;
	}
	public static final double[] clone(double[] in){
		int nx = in.length;
		double[]out = new double[nx];
		for(int i=0; i<nx; i++){
			out[i] = in[i];
		}
		return out;
	}
	public static final int[][] clone(int[][] in){
		int nx = in.length;
		int ny = in[0].length;
		int[][]out = new int[nx][ny];
		for(int i=0; i<nx; i++)for(int j=0; j<ny; j++){
			out[i][j] = in[i][j];
		}
		return out;
	}
	public static final float[][] clone(float[][] in){
		int nx = in.length;
		int ny = in[0].length;
		float[][]out = new float[nx][ny];
		for(int i=0; i<nx; i++)for(int j=0; j<ny; j++){
			out[i][j] = in[i][j];
		}
		return out;
	}
	public static final double[][] clone(double[][] in){
		int nx = in.length;
		int ny = in[0].length;
		double[][]out = new double[nx][ny];
		for(int i=0; i<nx; i++)for(int j=0; j<ny; j++){
			out[i][j] = in[i][j];
		}
		return out;
	}
	public static final int[][][] clone(int[][][] in){
		int nx = in.length;
		int ny = in[0].length;
		int nz = in[0][0].length;
		int[][][] out = new int[nx][ny][nz];
		for(int i=0; i<nx; i++)for(int j=0; j<ny; j++)for(int k=0; k<nz; k++){
			out[i][j][k] = in[i][j][k];
		}
		return out;
	}
	public static final double[][][] clone(double[][][] in){
		int nx = in.length;
		int ny = in[0].length;
		int nz = in[0][0].length;
		double[][][] out = new double[nx][ny][nz];
		for(int i=0; i<nx; i++)for(int j=0; j<ny; j++)for(int k=0; k<nz; k++){
			out[i][j][k] = in[i][j][k];
		}
		return out;
	}
	public static final float[][][] clone(float[][][] in){
		int nx = in.length;
		int ny = in[0].length;
		int nz = in[0][0].length;
		float[][][] out = new float[nx][ny][nz];
		for(int i=0; i<nx; i++)for(int j=0; j<ny; j++)for(int k=0; k<nz; k++){
			out[i][j][k] = in[i][j][k];
		}
		return out;
	}
	public static final double[][][][] clone(double[][][][] in){
		int nx = in.length;
		int ny = in[0].length;
		int nz = in[0][0].length;
		int nc = in[0][0][0].length;
		double[][][][] out = new double[nx][ny][nz][nc];
		for(int i=0; i<nx; i++)for(int j=0; j<ny; j++)for(int k=0; k<nz; k++)for(int l=0; l<nc; l++){
			out[i][j][k][l] = in[i][j][k][l];
		}
		return out;
	}
	public static final float[][][][] clone(float[][][][] in){
		int nx = in.length;
		int ny = in[0].length;
		int nz = in[0][0].length;
		int nc = in[0][0][0].length;
		float[][][][] out = new float[nx][ny][nz][nc];
		for(int i=0; i<nx; i++)for(int j=0; j<ny; j++)for(int k=0; k<nz; k++)for(int l=0; l<nc; l++){
			out[i][j][k][l] = in[i][j][k][l];
		}
		return out;
	}
	
	public static final double[][] toDouble(float[][] a){
		double[][] b = new double[a.length][a[0].length];
		for(int i=0; i<a.length; i++)for(int j=0; j<a[0].length; j++){
			b[i][j]=a[i][j];
		}
		return b;
	}
	public static final double[][] toDouble(int[][] a){
		double[][] b = new double[a.length][a[0].length];
		for(int i=0; i<a.length; i++)for(int j=0; j<a[0].length; j++){
			b[i][j]=a[i][j];
		}
		return b;
	}
	public static final float[][] toFloat(byte[][] a){
		float[][] b = new float[a.length][a[0].length];
		for(int i=0; i<a.length; i++)for(int j=0; j<a[0].length; j++){
			b[i][j]= (float)a[i][j];
		}
		return b;
	}
	public static final float[][][] toFloat(byte[][][] a){
		float[][][] b = new float[a.length][a[0].length][a[0][0].length];
		for(int i=0; i<a.length; i++)for(int j=0; j<a[0].length; j++)for(int k=0; k<a[0][0].length; k++){
			b[i][j][k]= (float)a[i][j][k];
		}
		return b;
	}
	public static final float[][] toFloat(double[][] a){
		float[][] b = new float[a.length][a[0].length];
		for(int i=0; i<a.length; i++)for(int j=0; j<a[0].length; j++){
			b[i][j]= (float)a[i][j];
		}
		return b;
	}
	public static final float[][] toFloat(int[][] a){
		float[][] b = new float[a.length][a[0].length];
		for(int i=0; i<a.length; i++)for(int j=0; j<a[0].length; j++){
			b[i][j]=a[i][j];
		}
		return b;
	}
	
	//***//
	
	public static final double[] toDouble(float[] a){
		double[] b = new double[a.length];
		for(int i=0; i<a.length; i++){
			b[i]=a[i];
		}
		return b;
	}
	public static final double[] toDouble(int[] a){
		double[] b = new double[a.length];
		for(int i=0; i<a.length; i++){
			b[i]=a[i];
		}
		return b;
	}
	
	public static final String[] toString(double[] a){
		String[] b = new String[a.length];
		for(int i=0; i<a.length; i++){
			b[i]=Double.toString(a[i]);
		}
		return b;
	}
	public static final String[] toString(float[] a){
		String[] b = new String[a.length];
		for(int i=0; i<a.length; i++){
			b[i]=Float.toString(a[i]);
		}
		return b;
	}
	public static final float[] toFloat(double[] a){
		float[] b = new float[a.length];
		for(int i=0; i<a.length; i++){
			b[i]=(float)a[i];
		}
		return b;
	}
	public static final int[] toInt(long[] a){
		int[] b = new int[a.length];
		for(int i=0; i<a.length; i++){
			b[i]= (int)a[i];
		}
		return b;
	}
	public static final long[] toLong(int[] a){
		long[] b = new long[a.length];
		for(int i=0; i<a.length; i++){
			b[i]=a[i];
		}
		return b;
	}
	public static final float[] toFloat(int[] a){
		float[] b = new float[a.length];
		for(int i=0; i<a.length; i++){
			b[i]=a[i];
		}
		return b;
	}
	public static final float[] toFloat(boolean[] a){
		float[] b = new float[a.length];
		for(int i=0; i<a.length; i++){
			
			if(a[i]) 
				b[i]=1f;
			
		}
		return b;
	}
	public static final int[] toInt(double[] a){
		int[] b = new int[a.length];
		for(int i=0; i<a.length; i++){
			b[i]=(int)a[i];
		}
		return b;
	}
	public static final int[] toInt(float[] a){
		int[] b = new int[a.length];
		for(int i=0; i<a.length; i++){
			b[i]=(int)a[i];
		}
		return b;
	}
	public static final int[] toIntRound(double[] a){
		int[] b = new int[a.length];
		for(int i=0; i<a.length; i++){
			b[i]=(int)Math.round(a[i]);
		}
		return b;
	}
	public static final int[] toIntRound(float[] a){
		int[] b = new int[a.length];
		for(int i=0; i<a.length; i++){
			b[i]=(int)Math.round(a[i]);
		}
		return b;
	}
	public static final int[][] toInt(double[][] a){
		int[][] b = new int[a.length][a[0].length];
		for(int i=0; i<a.length; i++)for(int j=0; j<a[0].length; j++){
			b[i][j]= (int)a[i][j];
		}
		return b;
	}
	public static final int[][] toInt(float[][] a){
		int[][] b = new int[a.length][a[0].length];
		for(int i=0; i<a.length; i++)for(int j=0; j<a[0].length; j++){
			b[i][j]= (int)a[i][j];
		}
		return b;
	}
	public static final int[][][] toInt(boolean[][][] a){
		int[][][] b = new int[a.length][a[0].length][a[0][0].length];
		for(int i=0; i<a.length; i++)for(int j=0; j<a[0].length; j++)for(int k=0; k<a[0][0].length; k++){
			if(a[i][j][k]) b[i][j][k]=1;
		}
		return b;
	}
	public static final float[] normalizeSum(float[] in){
		float sum = 0;
		float[] out = new float[in.length];
		for(int i=0; i<in.length; i++) {	sum+=in[i]; 		}
		for(int i=0; i<in.length; i++) {	out[i]=in[i]/sum; 	}
		return out;
	}
	public static final double[] normalizeSum(double[] in){
		double sum = 0;
		double[] out = new double[in.length];
		for(int i=0; i<in.length; i++) {	sum+=in[i]; 		}
		for(int i=0; i<in.length; i++) {	out[i]=in[i]/sum; 	}
		return out;
	}
	public static final float[] normalizeLength(float[] in){
		float[] out = new float[in.length];
		float sum = (float)Math.sqrt(sumSquares(in));
		for(int i=0; i<in.length; i++) {	out[i]=in[i]/sum; 	}
		return out;
	}
	public static final void normalizeComponentsOverwrite(float[][][][] in){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nc=in[0][0][0].length;
		float lensqr = -1;
		for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++){
			lensqr = sumSquares(in[x][y][z]);
			if(lensqr>0){
				for(int c=0; c<nc; c++){ in[x][y][z][c] = (float)(in[x][y][z][c]/Math.sqrt(lensqr)); }
			}
		}
	}
	static public float[][][] vectorMagnitude4D(float[][][][]  vec4D) {
		int rows=vec4D.length;
		int cols=vec4D[0].length;
		int slices=vec4D[0][0].length;
		int components = vec4D[0][0][0].length;
		float[][][] M = new float[rows][cols][slices];
		double sum = 0;
		for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) for (int k = 0; k < slices; k++) {
			sum = 0;
			for (int l = 0; l < components; l++) {
				sum += Math.pow(vec4D[i][j][k][l], 2);
			}
			M[i][j][k] = (float) Math.sqrt(sum);
		}
		return M;
	}
	public static final double[] normalizeLength(double[] in){
		double[] out = new double[in.length];
		double sum = Math.sqrt(sumSquares(in));
		for(int i=0; i<in.length; i++) {	out[i]=in[i]/sum; 	}
		return out;
	}
	public static final void normalizeLengthInPlace(double[] in){
		double sum = sumSquares(in);
		for(int i=0; i<in.length; i++) {	in[i] = in[i]/sum; 	}
	}
	
	public static final float[] normalizeTo(float[] in, float normsum){
		float sum = 0;
		float[] out = new float[in.length];
		for(int i=0; i<in.length; i++) {	sum+=in[i]; 		}
		for(int i=0; i<in.length; i++) {	out[i]=normsum*in[i]/sum; 	}
		return out;
	}
	public static final double[] normalizeTo(double[] in, double normsum){
		double sum = 0;
		double[] out = new double[in.length];
		for(int i=0; i<in.length; i++) {	sum+=in[i]; 		}
		for(int i=0; i<in.length; i++) {	out[i]=normsum*in[i]/sum; 	}
		return out;
	}
	public static final void normalizeToInPlace(double[] in, double normsum){
		double sum = 0;
		for(int i=0; i<in.length; i++) {	sum+=in[i]; 		}
		for(int i=0; i<in.length; i++) {	in[i] = normsum*in[i]/sum; 	}
	}
	
	public static final int max(int[] in){
		int out = in[0];
		for(int i=1; i<in.length; i++){
			if(in[i]>out){ out = in[i]; }
		}
		return out;
	}
	public static final double max(double[] in){
		double out = in[0];
		for(int i=1; i<in.length; i++){
			if(in[i]>out){ out = in[i]; }
		}
		return out;
	}
	public static final float max(float[] in){
		float out = in[0];
		for(int i=1; i<in.length; i++){
			if(in[i]>out){ out = in[i]; }
		}
		return out;
	}
	
	public static final int min(int[] in){
		int out = in[0];
		for(int i=1; i<in.length; i++){
			if(in[i]<out){ out = in[i]; }
		}
		return out;
	}
	public static final double min(double[] in){
		double out = in[0];
		for(int i=1; i<in.length; i++){
			if(in[i]<out){ out = in[i]; }
		}
		return out;
	}
	public static final float min(float[] in){
		float out = in[0];
		for(int i=1; i<in.length; i++){
			if(in[i]<out){ out = in[i]; }
		}
		return out;
	}
	
	public static final int maxIdx(int[] in){
		int out = in[0];
		for(int i=1; i<in.length; i++){
			if(in[i]>in[out]){ out = i; }
		}
		return out;
	}
	public static final int maxIdx(double[] in){
		int out = 0;
		for(int i=1; i<in.length; i++){
			if(in[i]>in[out]){ out = i; }
		}
		return out;
	}
	public static final int maxIdx(float[] in){
		int out = 0;
		for(int i=1; i<in.length; i++){
			if(in[i]>in[out]){ out = i; }
		}
		return out;
	}
	
	public static final int minIdx(int[] in){
		int out = in[0];
		for(int i=1; i<in.length; i++){
			if(in[i]<in[out]){ out = i; }
		}
		return out;
	}
	public static final int minIdx(double[] in){
		int out = 0;
		for(int i=1; i<in.length; i++){
			if(in[i]<in[out]){ out = i; }
		}
		return out;
	}
	public static final int minIdx(float[] in){
		int out = 0;
		for(int i=1; i<in.length; i++){
			if(in[i]<in[out]){ out = i; }
		}
		return out;
	}
	
	public static final int[] add(int[] a, int val){
		int[] out = new int[a.length];
		for(int i=0; i<a.length; i++) {out[i]=a[i]+val;}
		return out;
	}
	public static final long[] add(long[] a, long val){
		long[] out = new long[a.length];
		for(int i=0; i<a.length; i++) {out[i]=a[i]+val;}
		return out;
	}
	public static final float[] add(float[] a, float val){
		float[] out = new float[a.length];
		for(int i=0; i<a.length; i++) {out[i]=a[i]+val;}
		return out;
	}
	public static final double[] add(double[] a, double val){
		double[] out = new double[a.length];
		for(int i=0; i<a.length; i++) {out[i]=a[i]+val;}
		return out;
	}
	
	public static final void addInPlace(int[] a, int val){
		for(int i=0; i<a.length; i++) { a[i] += val; }
	}
	public static final void addInPlace(long[] a, long val){
		for(int i=0; i<a.length; i++) { a[i] += val; }
	}
	public static final void addInPlace(float[] a, float val){
		for(int i=0; i<a.length; i++) { a[i] += val; }
	}
	public static final void addInPlace(double[] a, double val){
		for(int i=0; i<a.length; i++) { a[i] += val; }
	}
	
	public static final int[] add(int[] a, int[] b){
		if(a.length!=b.length) return null;
		int[] out = new int[a.length];
		for(int i=0; i<b.length; i++) {out[i]=a[i]+b[i];}
		return out;
	}
	public static final float[] add(float[] a, float[] b){
		if(a.length!=b.length) return null;
		float[] out = new float[a.length];
		for(int i=0; i<b.length; i++) {out[i]=a[i]+b[i];}
		return out;
		
	}
	public static final void addInPlace(double[] a, double[] b){
		if(a.length!=b.length) return;
		
		for(int i=0; i<b.length; i++) a[i]+=b[i]; 
		
	}
	public static final void addInPlace(double[][] a, double[][] b){
		if(a.length!=b.length) return;
		if(a[0].length!=b[0].length) return;
		
		for(int i=0; i<b.length; i++) for(int j=0; j<b[0].length; j++){
			a[i][j]+=b[i][j]; 
		}
	}
	public static final void addInPlace(double[][][] a, double[][][] b){
		if(a.length!=b.length) return;
		if(a[0].length!=b[0].length) return;
		if(a[0][0].length!=b[0][0].length) return;
		
		for(int i=0; i<b.length; i++) for(int j=0; j<b[0].length; j++)for(int k=0; k<b[0][0].length; k++){
			a[i][j][k]+=b[i][j][k]; 
		}
	}
	public static final void addInPlace(double[][][][] a, double[][][][] b){
		if(a.length!=b.length) return;
		if(a[0].length!=b[0].length) return;
		if(a[0][0].length!=b[0][0].length) return;
		
		for(int i=0; i<b.length; i++) for(int j=0; j<b[0].length; j++)for(int k=0; k<b[0][0].length; k++) for(int l=0; l<b[0][0][0].length; l++){
			a[i][j][k][l]+=b[i][j][k][l]; 
		}
	}
	
	public static final void addInPlace(float[] a, float[] b){
		if(a.length!=b.length) return;
		
		for(int i=0; i<b.length; i++) a[i]+=b[i]; 
		
	}
	public static final void addInPlace(float[][] a, float[][] b){
		if(a.length!=b.length) return;
		if(a[0].length!=b[0].length) return;
		
		for(int i=0; i<b.length; i++) for(int j=0; j<b[0].length; j++){
			a[i][j]+=b[i][j]; 
		}
	}
	public static final void addInPlace(float[][][] a, float[][][] b){
		if(a.length!=b.length) return;
		if(a[0].length!=b[0].length) return;
		if(a[0][0].length!=b[0][0].length) return;
		
		for(int i=0; i<b.length; i++) for(int j=0; j<b[0].length; j++)for(int k=0; k<b[0][0].length; k++){
			a[i][j][k]+=b[i][j][k]; 
		}
	}
	public static final void addInPlace(float[][][][] a, float[][][][] b){
		if(a.length!=b.length) return;
		if(a[0].length!=b[0].length) return;
		if(a[0][0].length!=b[0][0].length) return;
		
		for(int i=0; i<b.length; i++) for(int j=0; j<b[0].length; j++)for(int k=0; k<b[0][0].length; k++) for(int l=0; l<b[0][0][0].length; l++){
			a[i][j][k][l]+=b[i][j][k][l]; 
		}
	}
	
	public static final void subtractInPlace(double[] a, double[] b){
		for(int i=0; i<a.length; i++){
			a[i]-=b[i];
		}
	}
	public static final int[] subtract(int[] a, int[] b){
		int[] c = new int[a.length];
		for(int i=0; i<a.length; i++){
			c[i] = a[i]-b[i];
		}
		return c;
	}
	public static final long[] subtract(long[] a, long[] b){
		long[] c = new long[a.length];
		for(int i=0; i<a.length; i++){
			c[i] = a[i]-b[i];
		}
		return c;
	}
	public static final float[] subtract(float[] a, float[] b){
		float[] c = new float[a.length];
		for(int i=0; i<a.length; i++){
			c[i] = a[i]-b[i];
		}
		return c;
	}
	public static final double[] subtract(double[] a, double[] b){
      double[] c = new double[a.length];
		for(int i=0; i<a.length; i++){
			c[i] = a[i]-b[i];
		}
		return c;
	}
	
	public static final double[] multiply(double[] a, double[] b){
		double[] out = new double[a.length];
		for(int i=0; i<a.length; i++){
			out[i]=a[i]*b[i]; 
		}
		return out;
	}
	public static final float[] multiply(float[] a, float[] b){
		float[] out = new float[a.length];
		for(int i=0; i<a.length; i++){
			out[i]=a[i]*b[i]; 
		}
		return out;
	}
	public static final int[] multiply(int[] a, int[] b){
		int[] out = new int[a.length];
		for(int i=0; i<a.length; i++){
			out[i]=a[i]*b[i]; 
		}
		return out;
	}
	public static final long[] multiply(long[] a, long[] b){
		long[] out = new long[a.length];
		for(int i=0; i<a.length; i++){
			out[i]=a[i]*b[i]; 
		}
		return out;
	}
	
	public static final double[] multiply(double[] in, double val){
		double[] out = new double[in.length];
		for(int i=0; i<in.length; i++){
			out[i]=in[i]*val; 
		}
		return out;
	}
	public static final void multiplyInPlace(double[] img, double val){
		for(int i=0; i<img.length; i++){
			img[i]*=val; 
		}
	}
	public static final float[] multiply(float[] in, float val){
		float[] out = new float[in.length];
		for(int i=0; i<in.length; i++){
			out[i]=in[i]*val; 
		}
		return out;
	}
	public static final void multiplyInPlace(float[] img, float val){
		for(int i=0; i<img.length; i++){
			img[i]*=val; 
		}
	}
	public static final void multiply(float[][][][] img, float val){
		for(int i=0; i<img.length; i++) for(int j=0; j<img[0].length; j++)for(int k=0; k<img[0][0].length; k++) for(int l=0; l<img[0][0][0].length; l++){
			img[i][j][k][l]*=val; 
		}
	}
	public static final void multiply(float[][][][] img, double val){
		for(int i=0; i<img.length; i++) for(int j=0; j<img[0].length; j++)for(int k=0; k<img[0][0].length; k++) for(int l=0; l<img[0][0][0].length; l++){
			img[i][j][k][l]*=val; 
		}
	}
	public static final void multiply(float[][][] img, float val){
		for(int i=0; i<img.length; i++) for(int j=0; j<img[0].length; j++)for(int k=0; k<img[0][0].length; k++){
			img[i][j][k]*=val; 
		}
	}
	public static final void multiply(float[][][] img, double val){
		for(int i=0; i<img.length; i++) for(int j=0; j<img[0].length; j++)for(int k=0; k<img[0][0].length; k++){
			img[i][j][k]*=val; 
		}
	}
	public static final void divide(float[][][][] img, float val){
		for(int i=0; i<img.length; i++) for(int j=0; j<img[0].length; j++)for(int k=0; k<img[0][0].length; k++) for(int l=0; l<img[0][0][0].length; l++){
			img[i][j][k][l]/=val; 
		}
	}
	public static final void divide(float[][] img, float val){
		for(int i=0; i<img.length; i++) for(int j=0; j<img[0].length; j++){
			img[i][j]/=val; 
		}
	}
	public static final void divide(float[][][][] img, double val){
		for(int i=0; i<img.length; i++) for(int j=0; j<img[0].length; j++)for(int k=0; k<img[0][0].length; k++) for(int l=0; l<img[0][0][0].length; l++){
			img[i][j][k][l]/=val; 
		}
	}
	public static final void divide(float[][][] img, float val){
		for(int i=0; i<img.length; i++) for(int j=0; j<img[0].length; j++)for(int k=0; k<img[0][0].length; k++){
			img[i][j][k]/=val; 
		}
	}
	public static final void divide(float[][][] img, double val){
		for(int i=0; i<img.length; i++) for(int j=0; j<img[0].length; j++)for(int k=0; k<img[0][0].length; k++){
			img[i][j][k]/=val; 
		}
	}
	public static final void divide(double[][] img, double val){
		for(int i=0; i<img.length; i++) for(int j=0; j<img[0].length; j++){
			img[i][j]/=val; 
		}
	}
	public static final void divide(double[] img, double val){
		for(int i=0; i<img.length; i++) {
			img[i]/=val; 
		}
	}
	public static final void divide(float[] img, float val){
		for(int i=0; i<img.length; i++) {
			img[i]/=val; 
		}
	}
	public static final void divide(int[] img, int val){
		for(int i=0; i<img.length; i++) {
			img[i]/=val; 
		}
	}
	public static final void negateInPlace(float[] in){
		for(int i=0; i<in.length; i++) {	in[i]*=-1; }
	}
	public static final float[] negate(float[] in){
		float[] out = new float[in.length];
		for(int i=0; i<in.length; i++) {	out[i]*=-in[i]; }
		return out;
	}
	public static final int sum(int[] in){
		int sum = 0;
		for(int i=0; i<in.length; i++) {	sum+=in[i]; 		}
		return sum;
	}
	public static final long sum(long[] in){
		long sum = 0;
		for(int i=0; i<in.length; i++) {	sum+=in[i]; 		}
		return sum;
	}
	public static final float sum(float[] in){
		float sum = 0;
		for(int i=0; i<in.length; i++) {	sum+=in[i]; 		}
		return sum;
	}
	public static final double sum(double[] in){
		double sum = 0;
		for(int i=0; i<in.length; i++) {	sum+=in[i]; 		}
		return sum;
	}
	public static long nnz(boolean[] in){
		long sum = 0;
		for(int i=0; i<in.length; i++) {	
			if(in[i]) { sum++; }
		}
		return sum;
	}
	public static long nnz(double[] in){
		long sum = 0;
		for(int i=0; i<in.length; i++) {	
			if(in[i]!=0) { sum++; }
		}
		return sum;
	}
	public static long nnz(float[] in){
		long sum = 0;
		for(int i=0; i<in.length; i++) {	
			if(in[i]!=0) { sum++; }
		}
		return sum;
	}
	public static long nnz(int[] in){
		long sum = 0;
		for(int i=0; i<in.length; i++) {	
			if(in[i]!=0) { sum++; }
		}
		return sum;
	}
	public static long nnz(long[] in){
		long sum = 0;
		for(int i=0; i<in.length; i++) {	
			if(in[i]!=0) { sum++; }
		}
		return sum;
	}
	public static final float prod(float[] in){
		float prod = 1f;
		for(int i=0; i<in.length; i++) {	prod*=in[i]; 		}
		return prod;
	}
	public static final double prod(double[] in){
		double prod = 1;
		for(int i=0; i<in.length; i++) {	prod*=in[i]; 		}
		return prod;
	}
	public static final float sumSquares(float[] in){
		float sum = 0;
		for(int i=0; i<in.length; i++) {	sum+=(in[i]*in[i]); 		}
		return sum;
	}
	public static final double sumSquares(double[] in){
		float sum = 0;
		for(int i=0; i<in.length; i++) {	sum+=(in[i]*in[i]); 		}
		return sum;
	}
	public static final float sumabs(float[] in){
		float sum = 0;
		for(int i=0; i<in.length; i++) {	sum+=Math.abs(in[i]); 		}
		return sum;
	}
	public static final double sumabs(double[] in){
		double sum = 0;
		for(int i=0; i<in.length; i++) {	sum+=Math.abs(in[i]); 		}
		return sum;
	}
	public static final void complement(boolean[][][] a){
		int nx=a.length;
		int ny=a[0].length;
		int nz=a[0][0].length;
		for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
			a[x][y][z] = !a[x][y][z];
		}
		
	}
	public static final String printArray(Object[] in){
		if(in==null) return "null";
		String out = "";	
		for(int i=0; i<in.length; i++){
			out += in[i].toString() +" ; ";
		}
		return out;
	}
	public static final String printArray(String[] in){
		if(in==null) return "null";
		String out = "";	
		for(int i=0; i<in.length; i++){
			out += in[i] +" ; ";
		}
		return out;
	}
	public static final String printArray(boolean[] in){
		if(in==null) return "null";
		String out = "";	
		for(int i=0; i<in.length; i++){
			if(in[i])
				out += "1 ; ";
			else
				out += "0 ; ";
		}
		return out;
	}
	
	public static final String printArray(byte[] in){
		if(in==null) return "null";
		String out = "";
		for(int i=0; i<in.length; i++){
			out += in[i] +" ";
		}
		return out;
	}
	
	public static final String printArray(float[] in){
		if(in==null) return "null";
		String out = "";	
		for(int i=0; i<in.length; i++){
			out += in[i] +" ";
		}
		return out;
	}

	public static final String printArray(double[] in){
		if(in==null) return "null";
		String out = "";
		for(int i=0; i<in.length; i++){
			out += in[i] +" ";
		}
		return out;
	}

	public static final String printArray(int[] in){
		if(in==null) return "null";
		String out = "";
		for(int i=0; i<in.length; i++){
			out += in[i] +" ";
		}
		return out;
	}
	public static final String printArray(long[] in){
		if(in==null) return "null";
		String out = "";
		for(int i=0; i<in.length; i++){
			out += in[i] +" ";
		}
		return out;
	}
	
	public static final String printArray(boolean[][] in){
		if(in==null) return "null";
		String out = "";
		for(int i=0; i<in.length; i++){
			for(int j=0; j<in[0].length; j++){
				if(in[i][j])
					out += "1  ";
				else
					out += "0  ";
			}
			out +="\n";
		}
		return out;
	}
	public static final String printArray(byte[][] in){
		if(in==null) return "null";
		String out = "";
		for(int i=0; i<in.length; i++){
			for(int j=0; j<in[0].length; j++){
				out += in[i][j] +" ";
			}
			out +="\n";
		}
		return out;
	}
	public static final String printArray(float[][] in){
		if(in==null) return "null";
		String out = "";
		for(int i=0; i<in.length; i++){
			for(int j=0; j<in[0].length; j++){
				out += in[i][j] +" ";
			}
			out +="\n";
		}
		return out;
	}
	
	public static final String printArray(int[][] in){
		if(in==null) return "null";
		String out = "";
		for(int i=0; i<in.length; i++){
			for(int j=0; j<in[0].length; j++){
				out += in[i][j] +" ";
			}
			out +="\n";
		}
		return out;
	}
	
	public static final String printArray(double[][] in){
		return printArray(in, " ", "\n");
	}
	public static final String printArray(double[][] in, String colsep, String rowsep){
		if(in==null) return "null";
		String out = "";
		for(int i=0; i<in.length; i++){
			for(int j=0; j<in[0].length; j++){
				out += in[i][j] + colsep;
			}
			out += rowsep;
		}
		return out;
	}
	
	public static final String printArray(float[][][] in){
		if(in==null) return "null";
		String out = "";
		for(int k=0; k<in[0][0].length; k++){
			
			for(int i=0; i<in.length; i++){
				for(int j=0; j<in[0].length; j++){
					out += in[i][j][k] +" ";
				}
				out +="\n";
			}
			out += "****************\n";
		}
		return out;
	}
	
	public static final String printArray(boolean[][][] in){
		if(in==null) return "null";
		String out = "";
		for(int k=0; k<in[0][0].length; k++){
			
			for(int i=0; i<in.length; i++){
				for(int j=0; j<in[0].length; j++){
					out += in[i][j][k] +" ";
				}
				out +="\n";
			}
			out += "****************\n";
		}
		return out;
	}
	public static final String printArray(int[][][] in){
		if(in==null) return "null";
		String out = "";
		for(int k=0; k<in[0][0].length; k++){
			
			for(int i=0; i<in.length; i++){
				for(int j=0; j<in[0].length; j++){
					out += in[i][j][k] +" ";
				}
				out +="\n";
			}
			out += "****************\n";
		}
		return out;
	}
	public static final String printArray(double[][][] in){
		if(in==null) return "null";
		String out = "";
		for(int k=0; k<in[0][0].length; k++){
			
			for(int i=0; i<in.length; i++){
				for(int j=0; j<in[0].length; j++){
					out += in[i][j][k] +" ";
				}
				out +="\n";
			}
			out += "****************\n";
		}
		return out;
	}
	public static final String printsubArray(byte[][][] in, int x, int y){
		if(in==null) return "null";
		String out = "";
		for(int i=0; i<in.length; i++){
			out += in[x][y][i] +" ";
		}
		return out;
	}
	
	public static final double[] getRow(double[][] in, int r){
		double[] out = new double[in[0].length];
		for(int i=0; i<in[0].length; i++){
			out[i]=in[r][i];
		}
		return out;
	}
	
	public static final double[] getColumn(double[][] in, int c){
		double[] out = new double[in.length];
		for(int i=0; i<in.length; i++){
			out[i]=in[i][c];
		}
		return out;
	}
	
	public static final float[] getRow(float[][] in, int r){
		float[] out = new float[in[0].length];
		for(int i=0; i<in[0].length; i++){
			out[i]=in[r][i];
		}
		return out;
	}
	
	public static final float[] getColumn(float[][] in, int c){
		float[] out = new float[in.length];
		for(int i=0; i<in.length; i++){
			out[i]=in[i][c];
		}
		return out;
	}
	
	public static final int[] getRow(int[][] in, int r){
		int[] out = new int[in[0].length];
		for(int i=0; i<in[0].length; i++){
			out[i]=in[r][i];
		}
		return out;
	}
	
	public static final int[] getColumn(int[][] in, int c){
		int[] out = new int[in.length];
		for(int i=0; i<in.length; i++){
			out[i]=in[i][c];
		}
		return out;
	}
	
	/**
	 * 
	 * @param in Input array
	 * @param x
	 * @param y
	 * @param dim The dimension to vary
	 * @return
	 */
	public static final byte[] get1dSubArray(byte[][][] in, int x, int y, int dim){
		byte[] out = null;
		int len = -1;
		if(dim==0){
			len = in.length;
			out = new byte[len];
			for(int i=0; i<len; i++){
				out[i]=in[i][x][y];
			}
		}else if(dim==1){
			len = in[0].length;
			out = new byte[len];
			for(int i=0; i<len; i++){
				out[i]=in[x][i][y];
			}
		}else if(dim==2){
			len = in[0][0].length;
			out = new byte[len];
			for(int i=0; i<len; i++){
				out[i]=in[x][y][i];
			}
		}else{
			System.out.println("Invalid dim");
		}
		return out;
	}
	
	/**
	 * 
	 * @param in Input array
	 * @param x
	 * @param y
	 * @param dim The dimension to vary
	 * @return
	 */
	public static final int[] get1dSubArray(int[][][] in, int x, int y, int dim){
		int[] out = null;
		int len = -1;
		if(dim==0){
			len = in.length;
			out = new int[len];
			for(int i=0; i<len; i++){
				out[i]=in[i][x][y];
			}
		}else if(dim==1){
			len = in[0].length;
			out = new int[len];
			for(int i=0; i<len; i++){
				out[i]=in[x][i][y];
			}
		}else if(dim==2){
			len = in[0][0].length;
			out = new int[len];
			for(int i=0; i<len; i++){
				out[i]=in[x][y][i];
			}
		}else{
			System.out.println("Invalid dim");
		}
		return out;
	}
	
	/**
	 * 
	 * @param in Input array
	 * @param x
	 * @param y
	 * @param dim The dimension to vary
	 * @return
	 */
	public static final float[] get1dSubArray(float[][][] in, int x, int y, int dim){
		float[] out = null;
		int len = -1;
		if(dim==0){
			len = in.length;
			out = new float[len];
			for(int i=0; i<len; i++){
				out[i]=in[i][x][y];
			}
		}else if(dim==1){
			len = in[0].length;
			out = new float[len];
			for(int i=0; i<len; i++){
				out[i]=in[x][i][y];
			}
		}else if(dim==2){
			len = in[0][0].length;
			out = new float[len];
			for(int i=0; i<len; i++){
				out[i]=in[x][y][i];
			}
		}else{
			System.out.println("Invalid dim");
		}
		return out;
	}
	
	/**
	 * 
	 * @param in
	 * @param x 
	 * @param dim The dimension to hold constant
	 * @return
	 */
	public static final byte[][] get2dSubArray(byte[][][] in, int x, int dim){
		byte[][] out = null;
		int dim1 = -1;
		int dim2 = -1;
		if(dim==0){
			dim1 = in[0].length;
			dim2 = in[0][0].length;
			out = new byte[dim1][dim2];
			for(int i=0; i<dim1; i++) for(int j=0; j<dim2; j++){
				out[i][j] = in[x][i][j];
			}
		}else if(dim==2){
			dim1 = in.length;
			dim2 = in[0][0].length;
			out = new byte[dim1][dim2];
			for(int i=0; i<dim1; i++) for(int j=0; j<dim2; j++){
				out[i][j] = in[i][x][j];
			}
		}else if(dim==3){
			dim1 = in.length;
			dim2 = in[0].length;
			out = new byte[dim1][dim2];
			for(int i=0; i<dim1; i++) for(int j=0; j<dim2; j++){
				out[i][j] = in[i][j][x];
			}
		}else{
			System.out.println("Invalid dim");
		}
		return out;
	}
	
	/**
	 * 
	 * @return
	 */
	public static final double[][] get1dSubArray(double[][] in, int[] x, boolean getRows){
		int N = x.length;
		int maxind = max(x);
		double[][] out = null;
		if(getRows){
			if(maxind>=in.length){
				System.out.println("INDICES REQUESTED OUT OF BOUNDS OF INPUT ROWS");
				return null;
			}
			out = new double[N][in[0].length];
		}else{
			if(N>=in[0].length){
				System.out.println("INDICES REQUESTED OUT OF BOUNDS OF INPUT COLUMNS");
				return null;
			}
			out = new double[in.length][N];
		}
		
		
		if(getRows){
			for(int i=0; i<N; i++)for(int j=0; j<in[0].length; j++){
				out[i][j] = in[x[i]][j];
			}
		}else{
			for(int i=0; i<in.length; i++)for(int j=0; j<N; j++){
				out[i][j] = in[j][x[i]];
			}
		}
		
		
		return out;
	}
	
	public static final int[][] get1dSubArray(int[][] in, int[] x, boolean getRows){
		int N = x.length;
		int maxind = max(x);
		int[][] out = null;
		if(getRows){
			if(maxind>=in.length){
				System.out.println("INDICES REQUESTED OUT OF BOUNDS OF INPUT ROWS");
				return null;
			}
			out = new int[N][in[0].length];
		}else{
			if(N>=in[0].length){
				System.out.println("INDICES REQUESTED OUT OF BOUNDS OF INPUT COLUMNS");
				return null;
			}
			out = new int[in.length][N];
		}
		
		
		if(getRows){
			for(int i=0; i<N; i++)for(int j=0; j<in[0].length; j++){
				out[i][j] = in[x[i]][j];
			}
		}else{
			for(int i=0; i<in.length; i++)for(int j=0; j<N; j++){
				out[i][j] = in[j][x[i]];
			}
		}
		
		
		return out;
	}
	
	
	public static final float[][] get1dSubArray(float[][] in, int[] x, boolean getRows){
		int N = x.length;
		int maxind = max(x);
		float[][] out = null;
		if(getRows){
			if(maxind>=in.length){
				System.out.println("INDICES REQUESTED OUT OF BOUNDS OF INPUT ROWS");
				return null;
			}
			out = new float[N][in[0].length];
		}else{
			if(N>=in[0].length){
				System.out.println("INDICES REQUESTED OUT OF BOUNDS OF INPUT COLUMNS");
				return null;
			}
			out = new float[in.length][N];
		}
		
		
		if(getRows){
			for(int i=0; i<N; i++)for(int j=0; j<in[0].length; j++){
				out[i][j] = in[x[i]][j];
			}
		}else{
			for(int i=0; i<in.length; i++)for(int j=0; j<N; j++){
				out[i][j] = in[j][x[i]];
			}
		}
		
		
		return out;
	}
	
	public static final int[] get1dSubArray(int[][] in, int fullidx, int startsub, int end, boolean fullRowSubCols){
		int[] out = new int[end-startsub];
		int k=0;
		
		if(fullRowSubCols){
			
			for(int i=startsub; i<end; i++){
				out[k]=in[fullidx][i];
				k++;
			}
			
		}else{
			for(int i=startsub; i<end; i++){
				out[k]=in[i][fullidx];
				k++;
			}
		}
		
		return out;
	}
	
	public static final float[] get1dSubArray(float[][] in, int fullidx, int startsub, int end, boolean fullRowSubCols){
		float[] out = new float[end-startsub];
		int k=0;
		
		if(fullRowSubCols){
			
			for(int i=startsub; i<end; i++){
				out[k]=in[fullidx][i];
				k++;
			}
			
		}else{
			for(int i=startsub; i<end; i++){
				out[k]=in[i][fullidx];
				k++;
			}
		}
		
		return out;
	}
	
	public static final void setComponent(float[][][][] dest, float[][][] src, int c, int nx, int ny, int nz){
		for(int x=0; x<nx; x++)for(int y=0; y<ny; y++)for(int z=0; z<nz; z++){
			dest[x][y][z][c]=src[x][y][z];
		}
	}
	public static final void setComponent(int[][][][] dest, int[][][] src, int c, int nx, int ny, int nz){
		for(int x=0; x<nx; x++)for(int y=0; y<ny; y++)for(int z=0; z<nz; z++){
			dest[x][y][z][c]=src[x][y][z];
		}
	}
	
	public static final double[][][] pad(double[][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		double[][][] out = new double[nx+padxlo+padxhi][ny+padylo+padyhi][nz+padzlo+padzhi];
		for(int x=0; x<nx; x++)for(int y=0; y<ny; y++)for(int z=0; z<nz; z++){
			out[x+padxlo][y+padylo][z+padzlo] = in[x][y][z];
		}
		return out;
	}
	public static final double[][][] pad(double[][][] in, int padx, int pady, int padz){
		return pad(in,padx,padx,pady,pady,padz,padz);
	}
	public static final float[][][] pad(float[][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		float[][][] out = new float[nx+padxlo+padxhi][ny+padylo+padyhi][nz+padzlo+padzhi];
		for(int x=0; x<nx; x++)for(int y=0; y<ny; y++)for(int z=0; z<nz; z++){
			out[x+padxlo][y+padylo][z+padzlo] = in[x][y][z];
		}
		return out;
	}
	public static final float[][][] pad(float[][][] in, int padx, int pady, int padz){
		return pad(in,padx,padx,pady,pady,padz,padz);
	}
	public static final int[][][] pad(int[][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int[][][] out = new int[nx+padxlo+padxhi][ny+padylo+padyhi][nz+padzlo+padzhi];
		for(int x=0; x<nx; x++)for(int y=0; y<ny; y++)for(int z=0; z<nz; z++){
			out[x+padxlo][y+padylo][z+padzlo] = in[x][y][z];
		}
		return out;
	}
	public static final int[][][] pad(int[][][] in, int padx, int pady, int padz){
		return pad(in,padx,padx,pady,pady,padz,padz);
	}
	public static final byte[][][] pad(byte[][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		byte[][][] out = new byte[nx+padxlo+padxhi][ny+padylo+padyhi][nz+padzlo+padzhi];
		for(int x=0; x<nx; x++)for(int y=0; y<ny; y++)for(int z=0; z<nz; z++){
			out[x+padxlo][y+padylo][z+padzlo] = in[x][y][z];
		}
		return out;
	}
	public static final byte[][][] pad(byte[][][] in, int padx, int pady, int padz){
		return pad(in,padx,padx,pady,pady,padz,padz);
	}
	public static final double[][] pad(double[][] in, int padxlo, int padxhi, int padylo, int padyhi){
		int nx=in.length;
		int ny=in[0].length;
		double[][] out = new double[nx+padxlo+padxhi][ny+padylo+padyhi];
		for(int x=0; x<nx; x++)for(int y=0; y<ny; y++){
			out[x+padxlo][y+padylo] = in[x][y];
		}
		return out;
	}
	public static final double[][] pad(double[][] in, int padx, int pady){
		return pad(in,padx,padx,pady,pady);
	}
	public static final float[][] pad(float[][] in, int padxlo, int padxhi, int padylo, int padyhi){
		int nx=in.length;
		int ny=in[0].length;
		float[][] out = new float[nx+padxlo+padxhi][ny+padylo+padyhi];
		for(int x=0; x<nx; x++)for(int y=0; y<ny; y++){
			out[x+padxlo][y+padylo] = in[x][y];
		}
		return out;
	}
	public static final float[][] pad(float[][] in, int padx, int pady){
		return pad(in,padx,padx,pady,pady);
	}
	public static final int[][] pad(int[][] in, int padxlo, int padxhi, int padylo, int padyhi){
		int nx=in.length;
		int ny=in[0].length;
		int[][] out = new int[nx+padxlo+padxhi][ny+padylo+padyhi];
		for(int x=0; x<nx; x++)for(int y=0; y<ny; y++){
			out[x+padxlo][y+padylo] = in[x][y];
		}
		return out;
	}
	public static final int[][] pad(int[][] in, int padx, int pady){
		return pad(in,padx,padx,pady,pady);
	}
	public static final byte[][] pad(byte[][] in, int padxlo, int padxhi, int padylo, int padyhi){
		int nx=in.length;
		int ny=in[0].length;
		byte[][] out = new byte[nx+padxlo+padxhi][ny+padylo+padyhi];
		for(int x=0; x<nx; x++)for(int y=0; y<ny; y++){
			out[x+padxlo][y+padylo] = in[x][y];
		}
		return out;
	}
	public static final byte[][] pad(byte[][] in, int padx, int pady){
		return pad(in,padx,padx,pady,pady);
	}
	
	//TODO: pad rep
	/**
	 * Pads an image with replicates of the first / last slice
	 */
	public static final int[][][] padRepeat(int[][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		int xx=-1, yy=-1, zz=-1; // coordinate in original image
		
		int[][][] out = new int[nx_new][ny_new][nz_new];
		for(int x=0; x<nx_new; x++)for(int y=0; y<ny_new; y++)for(int z=0; z<nz_new; z++){
			xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
			yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
			zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);
			
			out[x][y][z] = in[xx][yy][zz];
		}
		return out;
	}
	/**
	 * Pads an image with replicates of the first / last slice
	 */
	public static final int[][][][] padRepeat(int[][][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nc=in[0][0][0].length;
		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		int xx=-1, yy=-1, zz=-1; // coordinate in original image
		
		int[][][][] out = new int[nx_new][ny_new][nz_new][nc];
		for(int x=0; x<nx_new; x++)for(int y=0; y<ny_new; y++)for(int z=0; z<nz_new; z++)for(int c=0; c<nc; c++){
			xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
			yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
			zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);
			
			out[x][y][z][c] = in[xx][yy][zz][c];
		}
		return out;
	}
	public static final int[] padRepeatReshape1D(int[][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		int xx=-1, yy=-1, zz=-1; // coordinate in original image
		
		int[] out = new int[nx_new*ny_new*nz_new];
		int k = 0;
		if(rowwise){
			for(int x=0; x<nx_new; x++)for(int y=0; y<ny_new; y++)for(int z=0; z<nz_new; z++){
				xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
				yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
				zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);

				out[k] = in[xx][yy][zz];
				k++;
			}
		}else{
			for(int z=0; z<nz_new; z++)for(int y=0; y<ny_new; y++)for(int x=0; x<nx_new; x++){
				xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
				yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
				zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);

				out[k] = in[xx][yy][zz];
				k++;
			}
		}
		return out;
	}
	public static final int[][] padRepeatReshape2D(int[][][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nc=in[0][0][0].length;
		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		int xx=-1, yy=-1, zz=-1; // coordinate in original image
		
		int[][] out = new int[nx_new*ny_new*nz_new][nc];
		int k = 0;
		if(rowwise){
			for(int c=0; c<nc; c++){
				k = 0;
				for(int x=0; x<nx_new; x++)for(int y=0; y<ny_new; y++)for(int z=0; z<nz_new; z++){
					xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
					yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
					zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);

					out[k][c] = in[xx][yy][zz][c];
					k++;
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				k = 0;
				for(int z=0; z<nz_new; z++)for(int y=0; y<ny_new; y++)for(int x=0; x<nx_new; x++){
					xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
					yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
					zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);

					out[k][c] = in[xx][yy][zz][c];
					k++;
				}
			}
		}
		return out;
	}
	/**
	 * Pads an image with replicates of the first / last slice
	 */
	public static final float[][][] padRepeat(float[][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		int xx=-1, yy=-1, zz=-1; // coordinate in original image
		
		float[][][] out = new float[nx_new][ny_new][nz_new];
		for(int x=0; x<nx_new; x++)for(int y=0; y<ny_new; y++)for(int z=0; z<nz_new; z++){
			xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
			yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
			zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);
			
			out[x][y][z] = in[xx][yy][zz];
		}
		return out;
	}
	/**
	 * Pads an image with replicates of the first / last slice
	 */
	public static final float[][][][] padRepeat(float[][][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nc=in[0][0][0].length;
		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		int xx=-1, yy=-1, zz=-1; // coordinate in original image
		
		float[][][][] out = new float[nx_new][ny_new][nz_new][nc];
		for(int x=0; x<nx_new; x++)for(int y=0; y<ny_new; y++)for(int z=0; z<nz_new; z++)for(int c=0; c<nc; c++){
			xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
			yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
			zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);
			
			out[x][y][z][c] = in[xx][yy][zz][c];
		}
		return out;
	}
	public static final float[][] padRepeatReshape2D(float[][][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nc=in[0][0][0].length;
		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		int xx=-1, yy=-1, zz=-1; // coordinate in original image
		
//		System.out.println("nx: " + nx);
//		System.out.println("ny: " + ny);
//		System.out.println("nz: " + nz);
//		System.out.println("nc: " + nc);
		
		float[][] out = new float[nx_new*ny_new*nz_new][nc];
		int k = 0;
		if(rowwise){
			for(int c=0; c<nc; c++){
				k = 0;
				for(int x=0; x<nx_new; x++)for(int y=0; y<ny_new; y++)for(int z=0; z<nz_new; z++){
					xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
					yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
					zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);

					out[k][c] = in[xx][yy][zz][c];
					k++;
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				k = 0;
				for(int z=0; z<nz_new; z++)for(int y=0; y<ny_new; y++)for(int x=0; x<nx_new; x++){
					xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
					yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
					zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);

					out[k][c] = in[xx][yy][zz][c];
					k++;
				}
			}
		}
		
		return out;
	}
	public static final float[][] padRepeatReshape2Dcompwise(float[][][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nc=in[0][0][0].length;
		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		int xx=-1, yy=-1, zz=-1; // coordinate in original image
		
//		System.out.println("nx: " + nx);
//		System.out.println("ny: " + ny);
//		System.out.println("nz: " + nz);
//		System.out.println("nc: " + nc);
		
		float[][] out = new float[nc][nx_new*ny_new*nz_new];
		int k = 0;
		if(rowwise){
			for(int c=0; c<nc; c++){
				k = 0;
				for(int x=0; x<nx_new; x++)for(int y=0; y<ny_new; y++)for(int z=0; z<nz_new; z++){
					xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
					yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
					zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);

					out[c][k] = in[xx][yy][zz][c];
					k++;
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				k = 0;
				for(int z=0; z<nz_new; z++)for(int y=0; y<ny_new; y++)for(int x=0; x<nx_new; x++){
					xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
					yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
					zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);

					out[c][k] = in[xx][yy][zz][c];
					k++;
				}
			}
		}
		
		return out;
	}
	
	/**
	 * Pads an image with replicates of the first / last slice
	 */
	public static final double[][][] padRepeat(double[][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		int xx=-1, yy=-1, zz=-1; // coordinate in original image
		
		double[][][] out = new double[nx_new][ny_new][nz_new];
		for(int x=0; x<nx_new; x++)for(int y=0; y<ny_new; y++)for(int z=0; z<nz_new; z++){
			xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
			yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
			zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);
			
			out[x][y][z] = in[xx][yy][zz];
		}
		return out;
	}
	/**
	 * Pads an image with replicates of the first / last slice
	 */
	public static final double[][][][] padRepeat(double[][][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nc=in[0][0][0].length;
		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		int xx=-1, yy=-1, zz=-1; // coordinate in original image
		
		double[][][][] out = new double[nx_new][ny_new][nz_new][nc];
		for(int x=0; x<nx_new; x++)for(int y=0; y<ny_new; y++)for(int z=0; z<nz_new; z++)for(int c=0; c<nc; c++){
			xx = getReplicatePadCoordinate(x, nx, padxlo, padxhi);
			yy = getReplicatePadCoordinate(y, ny, padylo, padyhi);
			zz = getReplicatePadCoordinate(z, nz, padzlo, padzhi);
			
			out[x][y][z][c] = in[xx][yy][zz][c];
		}
		return out;
	}
	
	//TODO: pad mir
	/**
	 * Pads an image with mirroring of the first / last slices
	 */
	public static final int[][][] padMirror(int[][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		int xx=-1, yy=-1, zz=-1; // coordinate in original image
		
		int[][][] out = new int[nx_new][ny_new][nz_new];
		for(int x=0; x<nx_new; x++)for(int y=0; y<ny_new; y++)for(int z=0; z<nz_new; z++){
			xx = getMirrorPadCoordinate(x, nx, padxlo, padxhi);
			yy = getMirrorPadCoordinate(y, ny, padylo, padyhi);
			zz = getMirrorPadCoordinate(z, nz, padzlo, padzhi);
			
			out[x][y][z] = in[xx][yy][zz];
		}
		return out;
	}
	
	//TODO: pad mir
	/**
	 * Pads an image with mirroring of the first / last slices
	 */
	public static final double[][][] padMirror(double[][][] in, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		int xx=-1, yy=-1, zz=-1; // coordinate in original image
		
		double[][][] out = new double[nx_new][ny_new][nz_new];
		for(int x=0; x<nx_new; x++)for(int y=0; y<ny_new; y++)for(int z=0; z<nz_new; z++){
			xx = getMirrorPadCoordinate(x, nx, padxlo, padxhi);
			yy = getMirrorPadCoordinate(y, ny, padylo, padyhi);
			zz = getMirrorPadCoordinate(z, nz, padzlo, padzhi);
			
			out[x][y][z] = in[xx][yy][zz];
		}
		return out;
	}
	
	public static final int getReplicatePadCoordinate(int x, int nxOrig, int padxlo, int padxhi){
		int xx = -1;
		if(x<padxlo){ xx = 0;
		}else if( (x-padxlo) >= nxOrig){ xx = nxOrig - 1;
		}else{ xx=x-padxlo; }
		
		return xx;
	}
	
	public static final int getMirrorPadCoordinate(int x, int nxOrig, int padxlo, int padxhi){
		int xx = -1;
		if(x<padxlo){ xx = (padxlo - x -1);
		}else if( (x-2) >= nxOrig){ xx = 2*nxOrig -x + 1;
		}else{ xx=x-padxlo; }
		
		return xx;
	}
	
	public static int[] indexToCoordinates(int xyz, int nx, int ny, int nz){

		if(xyz>nx*ny*nz){
			return null;
		}else if(xyz<0){
			return null;
		}

		int[] coord = new int[3];

		coord[2] = xyz%nz;
		xyz -= coord[2];
		xyz/=nz;

		coord[1] = xyz%ny;
		xyz -= coord[1];
		xyz/=ny;

		coord[0] = xyz%nx;

		return coord;
	}

	public static int coordsToIndex(int x, int y, int z, int nx, int ny, int nz, boolean rowwise){
		if(rowwise)	return z + nz*y + nz*ny*x;
		else  		return x + nx*y + nx*ny*z;
	}
	
	public static final int[][][] cropReshape3D(int[] in, int nx, int ny, int nz, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi, boolean rowwise ){
		int[][][] out = new int[nx][ny][nz];

		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		if(rowwise){
			for(int x=0; x<nx; x++)for(int y=0; y<ny; y++)for(int z=0; z<nz; z++){
				out[x][y][z] = in[coordsToIndex(x+padxlo, y+padylo, z+padzlo, nx_new, ny_new, nz_new, rowwise)];
			}
		}else{
			for(int z=0; z<nz; z++)for(int y=0; y<ny; y++)for(int x=0; x<nx; x++){
				out[x][y][z] = in[coordsToIndex(x+padxlo, y+padylo, z+padzlo, nx_new, ny_new, nz_new, rowwise)];
			}
		}
		
		return out;
	}
	public static final float[][][] cropReshape3D(float[] in, int nx, int ny, int nz, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi, boolean rowwise ){
		float[][][] out = new float[nx][ny][nz];

		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		if(rowwise){
			for(int x=0; x<nx; x++)for(int y=0; y<ny; y++)for(int z=0; z<nz; z++){
				out[x][y][z] = in[coordsToIndex(x+padxlo, y+padylo, z+padzlo, nx_new, ny_new, nz_new, rowwise)];
			}
		}else{
			for(int z=0; z<nz; z++)for(int y=0; y<ny; y++)for(int x=0; x<nx; x++){
				out[x][y][z] = in[coordsToIndex(x+padxlo, y+padylo, z+padzlo, nx_new, ny_new, nz_new, rowwise)];
			}
		}
		
		return out;
	}
	
	public static final float[][][][] cropReshape4D(float[][] in, int nx, int ny, int nz, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi, boolean rowwise ){
		int nc = in.length;
		float[][][][] out = new float[nx][ny][nz][nc];

		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		if(rowwise){
			for(int c=0; c<nc; c++){
				for(int x=0; x<nx; x++)for(int y=0; y<ny; y++)for(int z=0; z<nz; z++){
					out[x][y][z][c] = in[c][coordsToIndex(x+padxlo, y+padylo, z+padzlo, nx_new, ny_new, nz_new, rowwise)];
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				for(int z=0; z<nz; z++)for(int y=0; y<ny; y++)for(int x=0; x<nx; x++){
					out[x][y][z][c] = in[c][coordsToIndex(x+padxlo, y+padylo, z+padzlo, nx_new, ny_new, nz_new, rowwise)];
				}
			}
		}
		
		return out;
	}
	
	public static final int[][][][] cropReshape4D(int[][] in, int nx, int ny, int nz, int padxlo, int padxhi, int padylo, int padyhi, int padzlo, int padzhi, boolean rowwise ){
		int nc = in.length;
		int[][][][] out = new int[nx][ny][nz][nc];

		int nx_new = nx+padxlo+padxhi;
		int ny_new = ny+padylo+padyhi;
		int nz_new = nz+padzlo+padzhi;
		
		if(rowwise){
			for(int c=0; c<nc; c++){
				for(int x=0; x<nx; x++)for(int y=0; y<ny; y++)for(int z=0; z<nz; z++){
					out[x][y][z][c] = in[c][coordsToIndex(x+padxlo, y+padylo, z+padzlo, nx_new, ny_new, nz_new, rowwise)];
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				for(int z=0; z<nz; z++)for(int y=0; y<ny; y++)for(int x=0; x<nx; x++){
					out[x][y][z][c] = in[c][coordsToIndex(x+padxlo, y+padylo, z+padzlo, nx_new, ny_new, nz_new, rowwise)];
				}
			}
		}
		
		return out;
	}
	
	public static final float[][] crop(float[][] in, int cropxlo, int cropxhi, int cropylo, int cropyhi){
		int nx=in.length;
		int ny=in[0].length;
		float[][] out = new float[nx-cropxlo-cropxhi][ny-cropxlo-cropxhi];
		for(int x=cropxlo; x<(nx-cropxhi); x++)for(int y=cropylo; y<(ny-cropyhi); y++){
			out[x-cropxlo][y-cropylo] = in[x][y];
		}
		return out;
	}
	public static final float[][] crop(float[][] in, int cropx, int cropy){
		return crop(in,cropx,cropx,cropy,cropy);
	}
	public static final int[][] crop(int[][] in, int cropxlo, int cropxhi, int cropylo, int cropyhi){
		int nx=in.length;
		int ny=in[0].length;
		int[][] out = new int[nx-cropxlo-cropxhi][ny-cropxlo-cropxhi];
		for(int x=cropxlo; x<(nx-cropxhi); x++)for(int y=cropylo; y<(ny-cropyhi); y++){
			out[x-cropxlo][y-cropylo] = in[x][y];
		}
		return out;
	}
	public static final int[][] crop(int[][] in, int cropx, int cropy){
		return crop(in,cropx,cropx,cropy,cropy);
	}
	public static final byte[][] crop(byte[][] in, int cropxlo, int cropxhi, int cropylo, int cropyhi){
		int nx=in.length;
		int ny=in[0].length;
		byte[][] out = new byte[nx-cropxlo-cropxhi][ny-cropxlo-cropxhi];
		for(int x=cropxlo; x<(nx-cropxhi); x++)for(int y=cropylo; y<(ny-cropyhi); y++){
			out[x-cropxlo][y-cropylo] = in[x][y];
		}
		return out;
	}
	public static final byte[][] crop(byte[][] in, int cropx, int cropy){
		return crop(in,cropx,cropx,cropy,cropy);
	}
	
	
	
	
	public static final float[] reshape1D(float[][][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		float[] out = new float[nx*ny*nz];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				out[i]=in[x][y][z];
				i++;
			}
		}else{
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[i]=in[x][y][z];
				i++;
			}
		}
		return out;
	}
	public static final boolean[] reshape1D(boolean[][][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		boolean[] out = new boolean[nx*ny*nz];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				out[i]=in[x][y][z];
				i++;
			}
		}else{
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[i]=in[x][y][z];
				i++;
			}
		}
		return out;
	}
	public static final byte[] reshape1D(byte[][][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		byte[] out = new byte[nx*ny*nz];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				out[i]=in[x][y][z];
				i++;
			}
		}else{
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[i]=in[x][y][z];
				i++;
			}
		}
		return out;
	}
	
	public static final int[] reshape1D(int[][][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int[] out = new int[nx*ny*nz];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				out[i]=in[x][y][z];
				i++;
			}
		}else{
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[i]=in[x][y][z];
				i++;
			}
		}
		return out;
	}
	
	public static final double[] reshape1D(double[][][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		double[] out = new double[nx*ny*nz];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				out[i]=in[x][y][z];
				i++;
			}
		}else{
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[i]=in[x][y][z];
				i++;
			}
		}
		return out;
	}
	
	public static final byte[] reshape1D(byte[][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		byte[] out = new byte[nx*ny];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++){
				out[i]=in[x][y];
				i++;
			}
		}else{
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[i]=in[x][y];
				i++;
			}
		}
		return out;
	}
	public static final int[] reshape1D(int[][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int[] out = new int[nx*ny];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++){
				out[i]=in[x][y];
				i++;
			}
		}else{
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[i]=in[x][y];
				i++;
			}
		}
		return out;
	}
	public static final float[] reshape1D(float[][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		float[] out = new float[nx*ny];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) {
				out[i]=in[x][y];
				i++;
			}
		}else{
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[i]=in[x][y];
				i++;
			}
		}
		return out;
	}
	public static final double[] reshape1D(double[][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		double[] out = new double[nx*ny];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) {
				out[i]=in[x][y];
				i++;
			}
		}else{
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[i]=in[x][y];
				i++;
			}
		}
		return out;
	}
	public static final double[][] reshape2D(double[] in, int nx, int ny, boolean rowwise){
		double[][] out = new double[nx][ny];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) {
				out[x][y]=in[i];
				i++;
			}
		}else{
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[x][y]=in[i];
				i++;
			}
		}
		return out;
	}
	public static final float[][] reshape2D(float[] in, int nx, int ny, boolean rowwise){
		float[][] out = new float[nx][ny];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) {
				out[x][y]=in[i];
				i++;
			}
		}else{
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[x][y]=in[i];
				i++;
			}
		}
		return out;
	}
	public static final int[][] reshape2D(int[] in, int nx, int ny, boolean rowwise){
		int[][] out = new int[nx][ny];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) {
				out[x][y]=in[i];
				i++;
			}
		}else{
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[x][y]=in[i];
				i++;
			}
		}
		return out;
	}
	public static final byte[][] reshape2D(byte[] in, int nx, int ny, boolean rowwise){
		byte[][] out = new byte[nx][ny];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) {
				out[x][y]=in[i];
				i++;
			}
		}else{
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[x][y]=in[i];
				i++;
			}
		}
		return out;
	}
	public static final boolean[][] reshape2D(boolean[] in, int nx, int ny, boolean rowwise){
		boolean[][] out = new boolean[nx][ny];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) {
				out[x][y]=in[i];
				i++;
			}
		}else{
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[x][y]=in[i];
				i++;
			}
		}
		return out;
	}
	public static final float[][] reshape2D(float[][][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		float[][] out = new float[nx*ny][nz];
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) {
				int i=0;
				for(int c=0; c<nz; c++){
					out[i][c]=in[x][y][c];
					i++;
				}
			}
		}else{
			int i =0;
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				for(int c=0; c<nz; c++){
					out[i][c]=in[x][y][c];
					i++;
				}
			}
		}
		return out;
	}
	public static final float[][] reshape2D(float[][][][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nc=in[0][0][0].length;
		float[][] out = new float[nx*ny*nz][nc];
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				int i=0;
				for(int c=0; c<nc; c++){
					out[i][c]=in[x][y][z][c];
					i++;
				}
				
			}
		}else{
			int i =0;
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				for(int c=0; c<nc; c++){
					out[i][c]=in[x][y][z][c];
					i++;
				}
				
			}
		}
		return out;
	}
	public static final int[][] reshape2Dcompwise(int[][][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int[][] out = new int[nz][nx*ny];
		int i=0;
		if(rowwise){
			i=0;
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++){
				for(int c=0; c<nz; c++){
					out[c][i]=in[x][y][c];
				}
				i++;
			}

		}else{
			i=0;
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				for(int c=0; c<nz; c++){
					out[c][i]=in[x][y][c];
				}
				i++;
			}
		}
		return out;
	}
	public static final float[][] reshape2Dcompwise(float[][][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		float[][] out = new float[nz][nx*ny];
		int i=0;
		if(rowwise){
			i=0;
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++){
				for(int c=0; c<nz; c++){
					out[c][i]=in[x][y][c];
				}
				i++;
			}

		}else{
			i=0;
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				for(int c=0; c<nz; c++){
					out[c][i]=in[x][y][c];
				}
				i++;
			}
		}
		return out;
	}
	public static final int[][] reshape2Dcompwise(int[][][][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nc=in[0][0][0].length;
		int[][] out = new int[nc][nx*ny*nz];
		int i=0;
		if(rowwise){
			i=0;
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				for(int c=0; c<nc; c++){
					out[c][i]=in[x][y][z][c];
				}
				i++;
			}

		}else{
			i=0;
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				for(int c=0; c<nc; c++){
					out[c][i]=in[x][y][z][c];
				}
				i++;
			}
		}
		return out;
	}
	public static final float[][] reshape2Dcompwise(float[][][][] in, boolean rowwise){
		int nx=in.length;
		int ny=in[0].length;
		int nz=in[0][0].length;
		int nc=in[0][0][0].length;
		float[][] out = new float[nc][nx*ny*nz];
		int i=0;
		if(rowwise){
			i=0;
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				for(int c=0; c<nc; c++){
					out[c][i]=in[x][y][z][c];
				}
				i++;
			}

		}else{
			i=0;
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				for(int c=0; c<nc; c++){
					out[c][i]=in[x][y][z][c];
				}
				i++;
			}
		}
		return out;
	}
	public static final int[][] reshape2DInt(byte[] in, int nx, int ny, boolean rowwise){
		int[][] out = new int[nx][ny];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) {
				out[x][y]=in[i];
				i++;
			}
		}else{
			for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[x][y]=in[i];
				i++;
			}
		}
		return out;
	}
	
	public static final float[][][] reshape3D(float[] in, int nx, int ny, int nz, boolean rowwise){
		float[][][] out = new float[nx][ny][nz];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				out[x][y][z]=in[i];
				i++;
			}
		}else{
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[x][y][z]=in[i];
				i++;
			}
		}
		return out;
	}
	public static final int[][][] reshape3D(int[] in, int nx, int ny, int nz, boolean rowwise){
		int[][][] out = new int[nx][ny][nz];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				out[x][y][z]=in[i];
				i++;
			}
		}else{
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[x][y][z]=in[i];
				i++;
			}
		}
		return out;
	}
	public static final double[][][] reshape3D(double[] in, int nx, int ny, int nz, boolean rowwise){
		double[][][] out = new double[nx][ny][nz];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				out[x][y][z]=in[i];
				i++;
			}
		}else{
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[x][y][z]=in[i];
				i++;
			}
		}
		return out;
	}
	public static final byte[][][] reshape3D(byte[] in, int nx, int ny, int nz, boolean rowwise){
		byte[][][] out = new byte[nx][ny][nz];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				out[x][y][z]=in[i];
				i++;
			}
		}else{
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[x][y][z]=in[i];
				i++;
			}
		}
		return out;
	}
	public static final int[][][] reshape3DInt(byte[] in, int nx, int ny, int nz, boolean rowwise){
		int[][][] out = new int[nx][ny][nz];
		int i=0;
		if(rowwise){
			for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
				out[x][y][z]=in[i];
				i++;
			}
		}else{
			for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
				out[x][y][z]=in[i];
				i++;
			}
		}
		return out;
	}
	public static final byte[][][] reshape3D(byte[][] in, int nx, int ny, int nc, boolean rowwise){
		byte[][][] out = new byte[nx][ny][nc];
		int i=0;
		if(rowwise){
			for(int c=0; c<nc; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++)  {
					out[x][y][c]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				i=0;
				for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][c]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	public static final int[][][] reshape3D(int[][] in, int nx, int ny, int nc, boolean rowwise){
		int[][][] out = new int[nx][ny][nc];
		int i=0;
		if(rowwise){
			for(int c=0; c<nc; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++)  {
					out[x][y][c]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				i=0;
				for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][c]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	public static final float[][][] reshape3D(float[][] in, int nx, int ny, int nc, boolean rowwise){
		float[][][] out = new float[nx][ny][nc];
		int i=0;
		if(rowwise){
			for(int c=0; c<nc; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++)  {
					out[x][y][c]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				i=0;
				for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][c]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	
	public static byte[][][] reshape3DSubset(byte[][] in, int nx, int ny, int ncmin, int ncmax, boolean rowwise){
		byte[][][] out = new byte[nx][ny][ncmax-ncmin];
		int i=0;
		if(rowwise){
			for(int c=ncmin; c<ncmax; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++)  {
					out[x][y][c-ncmin]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=ncmin; c<ncmax; c++){
				i=0;
				for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][c-ncmin]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	public static int[][][] reshape3DSubset(int[][] in, int nx, int ny, int ncmin, int ncmax, boolean rowwise){
		int[][][] out = new int[nx][ny][ncmax-ncmin];
		int i=0;
		if(rowwise){
			for(int c=ncmin; c<ncmax; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++)  {
					out[x][y][c-ncmin]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=ncmin; c<ncmax; c++){
				i=0;
				for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][c-ncmin]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	public static int[][][] reshape3DSubsetInt(byte[][] in, int nx, int ny, int ncmin, int ncmax, boolean rowwise){
		int[][][] out = new int[nx][ny][ncmax-ncmin];
		int i=0;
		if(rowwise){
			for(int c=ncmin; c<ncmax; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++)  {
					out[x][y][c-ncmin]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=ncmin; c<ncmax; c++){
				i=0;
				for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][c-ncmin]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	
	public static final double[][][][] reshape4D(double[][] in, int nx, int ny, int nz, int nc, boolean rowwise){
		double[][][][] out = new double[nx][ny][nz][nc];
		int i=0;
		if(rowwise){
			for(int c=0; c<nc; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
					out[x][y][z][c]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				i=0;
				for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][z][c]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	public static final float[][][][] reshape4D(float[][] in, int nx, int ny, int nz, int nc, boolean rowwise){
		float[][][][] out = new float[nx][ny][nz][nc];
		int i=0;
		if(rowwise){
			for(int c=0; c<nc; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
					out[x][y][z][c]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				i=0;
				for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][z][c]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	public static int[][][][] reshape4D(int[][] in, int nx, int ny, int nz, int nc, boolean rowwise){
		int[][][][] out = new int[nx][ny][nz][nc];
		int i=0;
		if(rowwise){
			for(int c=0; c<nc; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
					out[x][y][z][c]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				i=0;
				for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][z][c]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	
	public static int[][][][] reshape4DSubset(int[][] in, int nx, int ny, int nz, int ncmin, int ncmax, boolean rowwise){
		int[][][][] out = new int[nx][ny][nz][ncmax-ncmin];
		int i=0;
		if(rowwise){
			for(int c=ncmin; c<ncmax; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
					out[x][y][z][c-ncmin]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=ncmin; c<ncmax; c++){
				i=0;
				for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][z][c-ncmin]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	
	
	public static final byte[][][][] reshape4D(byte[][] in, int nx, int ny, int nz, int nc, boolean rowwise){
		byte[][][][] out = new byte[nx][ny][nz][nc];
		int i=0;
		if(rowwise){
			for(int c=0; c<nc; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
					out[x][y][z][c]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				i=0;
				for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][z][c]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	
	public static byte[][][][] reshape4DSubset(byte[][] in, int nx, int ny, int nz, int ncmin, int ncmax, boolean rowwise){
		byte[][][][] out = new byte[nx][ny][nz][ncmax-ncmin];
		int i=0;
		if(rowwise){
			for(int c=ncmin; c<ncmax; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
					out[x][y][z][c-ncmin]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=ncmin; c<ncmax; c++){
				i=0;
				for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][z][c-ncmin]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	
	public static final int[][][][] reshape4DInt(byte[][] in, int nx, int ny, int nz, int nc, boolean rowwise){
		int[][][][] out = new int[nx][ny][nz][nc];
		int i=0;
		if(rowwise){
			for(int c=0; c<nc; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
					out[x][y][z][c]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=0; c<nc; c++){
				i=0;
				for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][z][c]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	
	public static int[][][][] reshape4DSubsetInt(byte[][] in, int nx, int ny, int nz, int ncmin, int ncmax, boolean rowwise){
		int[][][][] out = new int[nx][ny][nz][ncmax-ncmin];
		int i=0;
		if(rowwise){
			for(int c=ncmin; c<ncmax; c++){
				i=0;
				for(int x=0; x<nx; x++)  for(int y=0; y<ny; y++) for(int z=0; z<nz; z++) {
					out[x][y][z][c-ncmin]=in[c][i];
					i++;
				}
			}
		}else{
			for(int c=ncmin; c<ncmax; c++){
				i=0;
				for(int z=0; z<nz; z++)  for(int y=0; y<ny; y++) for(int x=0; x<nx; x++){
					out[x][y][z][c-ncmin]=in[c][i];
					i++;
				}
			}
		}
		return out;
	}
	
	
	public static final double[][] transpose(double[][] in){
		int nr=in.length;
		int nc=in[0].length;
		double[][] out = new double[nc][nr];
		for(int i=0; i<nr; i++)  for(int j=0; j<nc; j++) {
			out[j][i] = in[i][j];
		}
		return out;
	}
	public static final float[][] transpose(float[][] in){
		int nr=in.length;
		int nc=in[0].length;
		float[][] out = new float[nc][nr];
		for(int i=0; i<nr; i++)  for(int j=0; j<nc; j++) {
			out[j][i] = in[i][j];
		}
		return out;
	}
	public static float[][] symMatrix(float[] v){
		float[][] out = null;
		int k=0;
		if(v.length==4){
			out = new float[2][2];
			for(int i=0; i<2; i++)for(int j=i; j<2; j++){
				out[i][j] = v[k];
				if(i!=j) out[j][i] = v[k];
			}
		}else if(v.length==6){
			out = new float[3][3];
			for(int i=0; i<3; i++)for(int j=i; j<3; j++){
				out[i][j] = v[k];
				if(i!=j) out[j][i] = v[k];
			}
		}
		return out;
	}
	public static float[][] symMatrix(float[][] in, boolean upper){

		float[][] out = new float[in.length][in[0].length];
		if(upper){
			for(int i=0; i<in.length; i++)for(int j=i; j<in[0].length; j++){
				out[i][j] = in[i][j];
				if(i!=j) out[j][i] = in[i][j]; 
			}
		}else{
			for(int j=0; j<in[0].length; j++) for(int i=j; i<in.length; i++){
				out[i][j] = in[i][j];
				if(i!=j) out[j][i] = in[i][j]; 
			}
		}

		return out;
	}
	public static double[][] symMatrix(double[] v){
		double[][] out = null;
		int k=0;
		if(v.length==4){
			out = new double[2][2];
			for(int i=0; i<2; i++)for(int j=i; j<2; j++){
				out[i][j] = v[k];
				if(i!=j){ out[j][i] = v[k]; }
				k++;
			}
		}else if(v.length==6){
			out = new double[3][3];
			for(int i=0; i<3; i++)for(int j=i; j<3; j++){
				out[i][j] = v[k];
				if(i!=j){ out[j][i] = v[k]; } 
				k++;
			}
		}
		return out;
	}
	public static double[][] symMatrixDouble(float[] v){
		double[][] out = null;
		int k=0;
		if(v.length==4){
			out = new double[2][2];
			for(int i=0; i<2; i++)for(int j=i; j<2; j++){
				out[i][j] = v[k];
				if(i!=j){ out[j][i] = v[k]; }
				k++;
			}
		}else if(v.length==6){
			out = new double[3][3];
			for(int i=0; i<3; i++)for(int j=i; j<3; j++){
				out[i][j] = v[k];
				if(i!=j){ out[j][i] = v[k]; } 
				k++;
			}
		}
		return out;
	}
	public static double[][] symMatrixDouble(float[][] in, boolean upper){

		double[][] out = new double[in.length][in[0].length];
		if(upper){
			for(int i=0; i<in.length; i++)for(int j=i; j<in[0].length; j++){
				out[i][j] = in[i][j];
				if(i!=j) out[j][i] = in[i][j]; 
			}
		}else{
			for(int j=0; j<in[0].length; j++) for(int i=j; i<in.length; i++){
				out[i][j] = in[i][j];
				if(i!=j) out[j][i] = in[i][j]; 
			}
		}

		return out;
	}
	
	
	
	
	public static final Object[] arrayFromObjectList(List<Object> l){
		Object[] out = new Object[l.size()];
		int i=0;
		for(Object n: l) { out[i]=n; i++; }
		return out;
	}
	
	public static final String[] arrayFromStringList(List<String> l){
		String[] out = new String[l.size()];
		int i=0;
		for(String n: l) { out[i]=n; i++; }
		return out;
	}
	public static final int[] arrayFromIntegerList(List<Integer> l){
		int[] out = new int[l.size()];
		int i=0;
		for(int n: l) { out[i]=n; i++; }
		return out;
	}
	public static final double[] arrayFromDoubleList(List<Double> l){
		double[] out = new double[l.size()];
		int i=0;
		for(double n: l) { out[i]=n; i++; }
		return out;
	}
	public static final float[] arrayFromFloatList(List<Float> l){
		float[] out = new float[l.size()];
		int i=0;
		for(float n: l) { out[i]=n; i++; }
		return out;
	}
	public static final float[][] arrayFromfloatarrayList(List<float[]> l){
		float[][] out = new float[l.size()][l.get(0).length];
		int i=0;
		for(float[] n: l) { out[i]=n; i++; }
		return out;
	}
	public static final int[][] arrayFromintarrayList(List<int[]> l){
		int[][] out = new int[l.size()][l.get(0).length];
		int i=0;
		for(int[] n: l) { out[i]=n; i++; }
		return out;
	}
	public static final byte[][] arrayFrombytearrayList(List<byte[]> l){
		byte[][] out = new byte[l.size()][l.get(0).length];
		int i=0;
		for(byte[] n: l) { out[i]=n; i++; }
		return out;
	}
	
	

	public static void main(String[] args){
		/*
		int[][][] test = new int[2][2][2];
		test[0][0][0]=0;
		test[0][0][1]=1;
		test[0][1][0]=2;
		test[0][1][1]=3;
		test[1][0][0]=4;
		test[1][0][1]=5;
		test[1][1][0]=6;
		test[1][1][1]=7;
		
		System.out.println(printArray(get1dSubArray(test,0,0,0)));
		System.out.println(printArray(get1dSubArray(test,0,0,0)));
		*/
		
		int[][] testa = new int[3][5];
		testa[0][0]=1;  testa[0][1]=2;  testa[0][2]=3;  testa[0][3]=4;  testa[0][4]=5;
		testa[1][0]=10;  testa[1][1]=20;  testa[1][2]=30;  testa[1][3]=40;  testa[1][4]=50;
		testa[2][0]=100;  testa[2][1]=200;  testa[2][2]=300;  testa[2][3]=400;  testa[2][4]=500;
		
		System.out.println(printArray(testa));
		System.out.println(printArray(get1dSubArray(testa, 2, 2, 5, true)));
	}
}
