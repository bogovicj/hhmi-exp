package net.imglib2.algorithms.crack;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import jitk.spline.KernelTransform;
import jitk.spline.ThinPlateR2LogRSplineKernelTransform;
import edu.jhu.ece.iacl.utility.ArrayUtil;
import mpicbg.models.TransformMesh;

public class SimCrack {

	int[] dims;
	int ndims;
	
	ThinPlateR2LogRSplineKernelTransform tps;
	
	int nCrackPts;
	double[][] crackCtrlPts;
	double[][] crackCtrlOffsets;
	
	private int edgePtNumSamples = 11;
	double eps 			 = 0.01;
	
	protected static Logger logger = LogManager.getLogger(SimCrack.class.getName());
	
	public SimCrack(){}
	
	public SimCrack(int[] dims )
	{
		setDims(dims);
	}
	
	public SimCrack(int[] dims, double[][] linePoints, double[][] offsets )
	{
		setDims(dims);
	}
	
	public void setDims(int[] dims){
		this.dims = dims;
		ndims = dims.length;
	}
	
	public int getEdgePtNumSamples() {
		return edgePtNumSamples;
	}

	public void setEdgePtNumSamples(int edgePtNumSamples) {
		this.edgePtNumSamples = edgePtNumSamples;
	}
	
	public void buildXfm3d()
	{

		int Nface = numFacePts3d( edgePtNumSamples );
		int N = Nface + ( 4 * nCrackPts );

		double[][] tpsCtrlPtsSrc = facePts3d( dims, edgePtNumSamples);
		double[][] tpsCtrlPtsTgt = facePts3d( dims, edgePtNumSamples);

		// set boundary conditions

		// set crack control points 
		int k = Nface;  // starting index for control points 
		for( int i = 0; i < nCrackPts; i++ )  
		{
			for( int d = 0; d < ndims; d++ ) 
			{
				// src
				tpsCtrlPtsSrc[d][k] = crackCtrlPts[d][i] + crackCtrlOffsets[d][i] + eps;     
				tpsCtrlPtsSrc[d][k] = crackCtrlPts[d][i] - crackCtrlOffsets[d][i] - eps;          
				// tgt
				tpsCtrlPtsSrc[d][k] = crackCtrlPts[d][i] + crackCtrlOffsets[d][i] + eps;     
				tpsCtrlPtsSrc[d][k] = crackCtrlPts[d][i] - crackCtrlOffsets[d][i] - eps;          

				// src
				tpsCtrlPtsTgt[d][k] = crackCtrlPts[d][i] + crackCtrlOffsets[d][i] + 2*eps;     
				tpsCtrlPtsTgt[d][k] = crackCtrlPts[d][i] - crackCtrlOffsets[d][i] - 2*eps;          
				// tgt
				tpsCtrlPtsTgt[d][k] = crackCtrlPts[d][i] + crackCtrlOffsets[d][i] + 2*eps;     
				tpsCtrlPtsTgt[d][k] = crackCtrlPts[d][i] - crackCtrlOffsets[d][i] - 2*eps;          
			}
		}

	}

   public static int numFacePts3d( int sPerDim )
   {

		int N = 
				2 * (sPerDim * sPerDim) +
				2 * (sPerDim * (sPerDim - 2)) +
				2 * (sPerDim - 2) * (sPerDim - 2);

      return N;
   }
	
	//public static double[][] facePoints( int[] dims, int samplesPerFace )
	//{
	//	int ndims = dims.length;
	//	int nfaces = 2*ndims;
	//	int k = ndims-1;
	//	int totSamples = samplesPerFace;
	//	while( k > 1){
	//		totSamples *= samplesPerFace;
	//		k--;
	//	}
	//	totSamples *= nfaces;
	//	System.out.println("nEdgeSamples : " + totSamples);
	//	
	//	double[][] res = new double[totSamples][ndims]; 
	//	
	//	int[] pos = new int[ndims];
	//	
	//	for (int d=0; d<ndims; d++){
	//		// vary d and fix the other dimensions
	//		
	//		
	//	}
	//	
	//	return res;
	//}
	
	/**
	 * 
	 * @param dims the dimensionality of the volume
	 * @param sPerDim number of samples in each dimension
     * @param dest the destination array
	 * @return
	 */
	public static double[][] facePts3d( int[] dims, int sPerDim  )
	{
		
		if( dims.length != 3 ){
			logger.error("facePts3d can only work with 3d volumes");
			return null;
		}
		int N = numFacePts3d( sPerDim ); 
		double[][] dest = new double[3][N];
		
		double delx = (dims[0]-1) / (double)(sPerDim-1);
		double dely = (dims[1]-1) / (double)(sPerDim-1);
		double delz = (dims[2]-1) / (double)(sPerDim-1);
		
		
		int x=0, y=0, z=0;
		double xx=0, yy=0, zz=0;
		
		int k = 0;
		
		// fix z
		for( x=0; x<sPerDim; x++) for( y=0; y<sPerDim; y++){
			xx = delx * x;
			yy = dely * y;
			zz = 0;
			
			// add point 
			insert(dest, k++, xx, yy, zz);
			
			zz = (sPerDim - 1)*delz;
			// add point
			insert(dest, k++, xx, yy, zz);

		}
	
		// fix y
		for( x=0; x<sPerDim; x++) for( z=1; z<(sPerDim-1); z++){
			xx = delx * x;
			yy = 0;
			zz = delz * z;
			
			insert(dest, k++, xx, yy, zz);
			
			yy = (sPerDim-1)* dely;
			insert(dest, k++, xx, yy, zz);
			
		}
		
		// fix x
		for( y=1; y<(sPerDim-1); y++) for( z=1; z<(sPerDim-1); z++){
			xx = 0;
			yy = dely * y;
			zz = delz * z;
			
			insert(dest, k++, xx, yy, zz);
			
			xx = (sPerDim-1) * delx;
			insert(dest, k++, xx, yy, zz);
		}
		
		return dest;
	}
	
	public static void insert( double[][] pts, int k, double x, double y, double z){
		pts[0][k] = x;
		pts[1][k] = y;
		pts[2][k] = z;
	}
	
	/**
	 * 
	 * @param ptList N x 3 array storing (x,y) positions of control points and width at each control point
	 */
	public void addParametrizedCrack(double[][] ptList){
		
	}
	
//	/**
//	 * 
//	 * @param fixedY
//	 * @param resX
//	 * @param widthSlope
//	 * @return
//	 */
//	public double[][] genLinearCrackParams(float fixedY, int resX, float widthSlope){
//		double[][] crackParam = new double[resX][3];
//		
//		//float w = xfmMesh.getWidth();
//		//float h = xfmMesh.getWidth();
//		
//		double step = w/(resX-1);
//		for(int i=0; i<resX; i++){
//			crackParam[i][0] = i*step;
//			crackParam[i][1] = fixedY;
//			crackParam[i][2] = widthSlope * (i*step);
//		}
//		
//		return crackParam;
//	}
	
//	public static void firstPass(){
//		float width = 20f, height = 20f;
//		int meshResX = 4, meshResY = 4;
//		SimCrack sc = new SimCrack(meshResX, meshResY, width, height);
//		
//		double[][] p = sc.genLinearCrackParams(5.5f, 5, 0.5f);
//		
//		System.out.println("\n" + ArrayUtil.printArray(p));
//		
//	}
	
	/**
	 * Generates a crack in the x direction
	 * @param nx image width in x
	 * @param ny image width in y
	 * @param sx starting
	 * @param yx
	 * @param N
	 * @param crackWidth
	 */
	public void genSimpleCrackXZ3d( int nx, int nz, int sx, int crackLength, int fixedY, double crackWidth )
	{
		crackCtrlPts 	 = new double[3][crackLength];
		crackCtrlOffsets = new double[3][crackLength];

		// first point should have zero offset 

		for ( int z=0; z<nz; z++){
			for ( int x=sx; x<(sx+crackLength); x++){
				for ( int i=0; i<crackLength; i++){

					// offsets for first and last crack points should be zero
					// otherwise, set non-zero crack offsets
					if(i>0 && i<(crackLength-1)){
						crackCtrlOffsets[1][i] = crackWidth;  // only offset in y 
					}

					crackCtrlPts[0][i] = x;
					crackCtrlPts[1][i] = fixedY;
					crackCtrlPts[2][i] = z;
				}
			}
		}
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("starting\n");

		int[] dims = new int[]{10,10,5};
//		SimCrack sc = new SimCrack(dims);
//		sc.buildXfm();
		
		double[][] ptsOut = facePts3d( dims, 7 );
		
		System.out.println("\n" + ArrayUtil.printArray(
				ArrayUtil.transpose(ptsOut),
            " ", ";" 
				));
		
		
		
		
		System.out.println("\nfinished");
	}



}
