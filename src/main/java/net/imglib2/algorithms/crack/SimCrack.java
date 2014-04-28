package net.imglib2.algorithms.crack;

import ij.IJ;
import ij.ImagePlus;
import mpicbg.ij.ThinPlateSplineMapping;
import net.imglib2.exception.ImgLibException;
import net.imglib2.img.Img;
import net.imglib2.img.imageplus.ImagePlusImg;
import net.imglib2.img.imageplus.ImagePlusImgs;
import net.imglib2.type.numeric.integer.ByteType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import jitk.spline.KernelTransform;
import jitk.spline.ThinPlateR2LogRSplineKernelTransform;
import jitk.spline.ThinPlateR2LogRSplineKernelTransformFloatSep;
import edu.jhu.ece.iacl.utility.ArrayUtil;

//import mpicbg.models.TransformMesh;

public class SimCrack {

	int[] dims;
	int ndims;
	
	ThinPlateR2LogRSplineKernelTransformFloatSep tps;
	
	int nCrackPts;
	double[][] crackCtrlPts;
	double[][] crackCtrlOffsets;
	
	private int edgePtNumSamples = 21;
	double eps 			 = 0.1;
	
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

		System.out.println("nface " + Nface);
		System.out.println("nCrackPts " + nCrackPts);
		System.out.println("N " + N);
		
		System.out.println("ctrlpt sz " + crackCtrlPts.length + " "+ crackCtrlPts[0].length);
		System.out.println("ctrlpt sz " + crackCtrlOffsets.length + " "+ crackCtrlOffsets[0].length);
		
		
		float[][] tpsCtrlPtsSrc = new float[3][N];
		float[][] tpsCtrlPtsTgt = new float[3][N];
				
		facePts3dFloat( dims, edgePtNumSamples, tpsCtrlPtsSrc, 0 );
		facePts3dFloat( dims, edgePtNumSamples, tpsCtrlPtsTgt, 0 );

		// set boundary conditions

		// set crack control points 
		int k = Nface;  // starting index for control points 
		for( int i = 0; i < nCrackPts; i++ )  
		{
			for( int d = 0; d < ndims; d++ ) 
			{
				// +eps
				tpsCtrlPtsSrc[d][k] = (float)(crackCtrlPts[d][i] + eps);
				tpsCtrlPtsTgt[d][k] = (float)(crackCtrlPts[d][i] + crackCtrlOffsets[d][i] + eps);   
				
				// -eps
				tpsCtrlPtsSrc[d][k+1] = (float)(crackCtrlPts[d][i] - eps);
				tpsCtrlPtsTgt[d][k+1] = (float)(crackCtrlPts[d][i] - crackCtrlOffsets[d][i] - eps);   
				    
				// +2eps
				tpsCtrlPtsSrc[d][k+2] = (float)(crackCtrlPts[d][i] + 2*eps);     
				tpsCtrlPtsTgt[d][k+2] = (float)(crackCtrlPts[d][i] + crackCtrlOffsets[d][i] + 2*eps);          
				
				// -2eps
				tpsCtrlPtsSrc[d][k+3] = (float)(crackCtrlPts[d][i] - 2*eps);     
				tpsCtrlPtsTgt[d][k+3] = (float)(crackCtrlPts[d][i] - crackCtrlOffsets[d][i] - 2*eps);      
				
			}
			k+=4;
			
		}
		
		tps = new ThinPlateR2LogRSplineKernelTransformFloatSep(3);
		tps.setLandmarks( tpsCtrlPtsSrc, tpsCtrlPtsTgt );
//		tps.setLandmarks( tpsCtrlPtsTgt, tpsCtrlPtsSrc);
		tps.fit();
		
		for( int i = Nface; i < N; i++ )  
		{
			System.out.println("\nmatch " + i + ": \n" +
					tpsCtrlPtsSrc[0][i] + " " + tpsCtrlPtsSrc[1][i] + " " + tpsCtrlPtsSrc[2][i]
					+ "\n -> \n" +
					tpsCtrlPtsTgt[0][i] + " " + tpsCtrlPtsTgt[1][i] + " " + tpsCtrlPtsTgt[2][i]
			+"\n");
		}

	}
	public void buildXfm2d()
	{

		int Nface = numFacePts2d( edgePtNumSamples );
		int N = Nface + ( 4 * nCrackPts );

		System.out.println("nface " + Nface);
		System.out.println("nCrackPts " + nCrackPts);
		System.out.println("N " + N);
		
		System.out.println("ctrlpt sz " + crackCtrlPts.length + " "+ crackCtrlPts[0].length);
		System.out.println("ctrlpt sz " + crackCtrlOffsets.length + " "+ crackCtrlOffsets[0].length);
		
		
		float[][] tpsCtrlPtsSrc = new float[ndims][N];
		float[][] tpsCtrlPtsTgt = new float[ndims][N];
				
		facePts2dFloat( dims, edgePtNumSamples, tpsCtrlPtsSrc, 0 );
		facePts2dFloat( dims, edgePtNumSamples, tpsCtrlPtsTgt, 0 );

		// set boundary conditions

		// set crack control points 
		int k = Nface;  // starting index for control points 
		for( int i = 0; i < nCrackPts; i++ )  
		{
			for( int d = 0; d < ndims; d++ ) 
			{
				// +eps
				tpsCtrlPtsSrc[d][k] = (float)(crackCtrlPts[d][i] + eps);
				tpsCtrlPtsTgt[d][k] = (float)(crackCtrlPts[d][i] + crackCtrlOffsets[d][i] + eps);   
				
				// -eps
				tpsCtrlPtsSrc[d][k+1] = (float)(crackCtrlPts[d][i] - eps);
				tpsCtrlPtsTgt[d][k+1] = (float)(crackCtrlPts[d][i] - crackCtrlOffsets[d][i] - eps);   
				    
				// +2eps
				tpsCtrlPtsSrc[d][k+2] = (float)(crackCtrlPts[d][i] + crackCtrlOffsets[d][i] + 2*eps);     
				tpsCtrlPtsTgt[d][k+2] = (float)(crackCtrlPts[d][i] + crackCtrlOffsets[d][i] + 4*eps);          
				
				// -2eps
				tpsCtrlPtsSrc[d][k+3] = (float)(crackCtrlPts[d][i] - crackCtrlOffsets[d][i] - 2*eps);     
				tpsCtrlPtsTgt[d][k+3] = (float)(crackCtrlPts[d][i] - crackCtrlOffsets[d][i] - 4*eps);      
				
			}
			k+=4;
			
		}
		
		tps = new ThinPlateR2LogRSplineKernelTransformFloatSep(ndims);
//		tps.setLandmarks( tpsCtrlPtsSrc, tpsCtrlPtsTgt );
		tps.setLandmarks( tpsCtrlPtsTgt, tpsCtrlPtsSrc);
		tps.fit();
		
		for( int i = Nface; i < N; i++ )  
		{
			System.out.println("\nmatch " + i + ": \n" +
					tpsCtrlPtsSrc[0][i] + " " + tpsCtrlPtsSrc[1][i] 
					+ "\n -> \n" +
					tpsCtrlPtsTgt[0][i] + " " + tpsCtrlPtsTgt[1][i]
			+"\n");
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
	
   public static int numFacePts2d( int sPerDim )
   {

		int N = 
				2 * (sPerDim ) +
				2 * (sPerDim - 2) ;

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
	 */
	public static void facePts3d( int[] dims, int sPerDim, double[][] dest, int startingIdx  )
	{
		
		if( dims.length != 3 ){
			logger.error("facePts3d can only work with 3d volumes");
			return;
		}
		int N = numFacePts3d( sPerDim ); 
		
		double delx = (dims[0]-1) / (double)(sPerDim-1);
		double dely = (dims[1]-1) / (double)(sPerDim-1);
		double delz = (dims[2]-1) / (double)(sPerDim-1);
		
		
		int x=0, y=0, z=0;
		double xx=0, yy=0, zz=0;
		
		int k = startingIdx;
		
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
		
	}
	
	/**
	 * 
	 * @param dims the dimensionality of the volume
	 * @param sPerDim number of samples in each dimension
     * @param dest the destination array
	 */
	public static void facePts3dFloat( int[] dims, int sPerDim, float[][] dest, int startingIdx  )
	{
		
		if( dims.length != 3 ){
			logger.error("facePts3d can only work with 3d volumes");
			return;
		}
		int N = numFacePts3d( sPerDim ); 
		if(dest.length != 3){
			logger.error(" input dest matrix must be length 3 in first dimension");
		}
		if(dest[0].length < N){
			logger.error(" input dest matrix is too short, must be of length " + N + 
					" or greater in 2nd dimension.");
		}
		
		float delx = (dims[0]-1) / (float)(sPerDim-1);
		float dely = (dims[1]-1) / (float)(sPerDim-1);
		float delz = (dims[2]-1) / (float)(sPerDim-1);
		
		
		int x=0, y=0, z=0;
		float xx=0, yy=0, zz=0;
		
		int k = startingIdx;
		
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
		
	}
	
	
	/**
	 * 
	 * @param dims the dimensionality of the volume
	 * @param sPerDim number of samples in each dimension
     * @param dest the destination array
	 */
	public static void facePts2dFloat( int[] dims, int sPerDim, float[][] dest, int startingIdx  )
	{
		
		if( dims.length != 2 ){
			logger.error("facePts3d can only work with 2d volumes");
			return;
		}
		int N = numFacePts2d( sPerDim ); 
		if(dest.length != 2){
			logger.error(" input dest matrix must be length 2 in first dimension");
		}
		if(dest[0].length < N){
			logger.error(" input dest matrix is too short, must be of length " + N + 
					" or greater in 2nd dimension.");
		}
		
		float delx = (dims[0]-1) / (float)(sPerDim-1);
		float dely = (dims[1]-1) / (float)(sPerDim-1);
		
		
		int    x=0,  y=0;
		float xx=0, yy=0;
		
		int k = startingIdx;
		
	
		// fix y
		for( x=0; x<sPerDim; x++){
			xx = delx * x;
			yy = 0;
			
			insert(dest, k++, xx, yy ); // y = 0
			
			yy = (sPerDim-1)* dely;
			insert(dest, k++, xx, yy ); // y = 'maxy'
			
		}
		
		// fix x
		for( y=1; y<(sPerDim-1); y++) {
			xx = 0;
			yy = dely * y;
		
			insert(dest, k++, xx, yy ); // x = 0
			
			xx = (sPerDim-1) * delx;
			insert(dest, k++, xx, yy ); // x = 'maxx'
		}
		
	}
	
	/**
	 * Unsafe
	 */
	public static void insert( double[][] pts, int k, double x, double y, double z){
		pts[0][k] = x;
		pts[1][k] = y;
		pts[2][k] = z;
	}
	
	/**
	 * Unsafe
	 */
	public static void insert( float[][] pts, int k, float x, float y, float z){
		pts[0][k] = x;
		pts[1][k] = y;
		pts[2][k] = z;
	}
	
	/**
	 * Unsafe
	 */
	public static void insert( double[][] pts, int k, double x, double y ){
		pts[0][k] = x;
		pts[1][k] = y;
	}
	
	/**
	 * Unsafe
	 */
	public static void insert( float[][] pts, int k, float x, float y ){
		pts[0][k] = x;
		pts[1][k] = y;
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
	 *  Generates a crack in the x direction 
	 * @param nx image width in x
	 * @param startX starting x 
	 * @param crackLength
	 * @param fixedY
	 * @param crackHalfWidth 
	 */
	public void genSimpleCrackX2d( int nx, int startX, int crackLength, int step, int fixedY, double crackHalfWidth )
	{
		crackCtrlPts 	 = new double[2][crackLength];
		crackCtrlOffsets = new double[2][crackLength];
		nCrackPts = 0;

		// first point should have zero offset 

			for ( int i=0; i<crackLength; i+=step){

				// offsets for first and last crack points should be zero
				// otherwise, set non-zero crack offsets
				if(i>0 && i<(crackLength-1)){
					crackCtrlOffsets[1][nCrackPts] = crackHalfWidth;  // only offset in y 
				}

				crackCtrlPts[0][nCrackPts] = startX + i;
				crackCtrlPts[1][nCrackPts] = fixedY;
				
				nCrackPts++;
			}
		
		System.out.println("crack pts: \n" + ArrayUtil.printArray(crackCtrlPts));
		System.out.println("\ncrack offs: \n" + ArrayUtil.printArray(crackCtrlOffsets));
	}
	
	
	/**
	 * Generates a crack in the x direction
	 * @param nx image width in x
	 * @param nz image width in z
	 * @param startX starting
	 * @param yx
	 * @param N
	 * @param crackWidth
	 */
	public void genSimpleCrackXZ3d( int nx, int nz, int startX, int crackLength, int fixedY, double crackHalfWidth )
	{
		crackCtrlPts 	 = new double[3][crackLength*nz];
		crackCtrlOffsets = new double[3][crackLength*nz];
		nCrackPts = 0;

		// first point should have zero offset 

		for ( int z=0; z<nz; z++){
			for ( int i=0; i<crackLength; i++){

				
				// offsets for first and last crack points should be zero
				// otherwise, set non-zero crack offsets
				if(i>0 && i<(crackLength-1)){
					crackCtrlOffsets[1][nCrackPts] = crackHalfWidth;  // only offset in y 
				}

				crackCtrlPts[0][nCrackPts] = startX + i;
				crackCtrlPts[1][nCrackPts] = fixedY;
				crackCtrlPts[2][nCrackPts] = z;
				
				nCrackPts++;
			}

		}
		
		System.out.println("crack pts: \n" + ArrayUtil.printArray(crackCtrlPts));
		System.out.println("\ncrack offs: \n" + ArrayUtil.printArray(crackCtrlOffsets));
	}
	
	public static void test3d() {
		int[] dims = new int[]{15,15,3};
		SimCrack sc = new SimCrack();
		sc.setDims(dims);
		
		int startX 		= 5;
		int crackLength = 5;
		int fixedY 		= 8;
		
		double crackHalfWidth = 2.0;
		sc.genSimpleCrackXZ3d(dims[0], dims[2], startX, crackLength, fixedY, crackHalfWidth);
		sc.buildXfm3d();
		
		System.out.println("sc: " + sc);
		
		FloatType b = new FloatType();
		Img<FloatType> img = ImgOps.createGradientImgY(dims[0], dims[1], dims[2], b);
		
		ImagePlus ip = null;
		try {
			ip = ImgOps.copyToImagePlus(img).getImagePlus();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		ImagePlus dest = new ImagePlus("dest", ip.getStack());
		
		System.out.println("in size:  " + ip.getWidth() + " " +  ip.getHeight() + " " + ip.getStackSize());
		System.out.println("dest size:  " + dest.getWidth() + " " +  dest.getHeight() + " " + dest.getStackSize());
		
		
		float[] testpt1 = new float[]{7, 7.0f, 1};
		float[] testpt2 = new float[]{7, 7.9f, 1};
		float[] testpt3 = new float[]{7, 8.1f, 1};
		float[] testpt4 = new float[]{7, 9.0f, 1};
		
		
		
		float[] outpt1 = sc.tps.transformPoint(testpt1);
		float[] outpt2 = sc.tps.transformPoint(testpt2);
		float[] outpt3 = sc.tps.transformPoint(testpt3);
		float[] outpt4 = sc.tps.transformPoint(testpt4);
		
		System.out.println("outpt1:\n" + ArrayUtil.printArray(outpt1));
		System.out.println("outpt2:\n" + ArrayUtil.printArray(outpt2));
		System.out.println("outpt3:\n" + ArrayUtil.printArray(outpt3));
		System.out.println("outpt4:\n" + ArrayUtil.printArray(outpt4));
		
		
		
		ThinPlateSplineMapping.mapInterval(sc.tps, ip, dest);
		
		String basename = "/groups/jain/home/bogovicj/projects/crackPatching/toyData/";
		
		IJ.save(ip, basename + "gradyImg.tif");
		IJ.save(dest, basename + "gradyImg_crack.tif");
		
		
//		double[][] ptsOut = facePts3d( dims, 7 );
//		
//		System.out.println("\n" + ArrayUtil.printArray(
//				ArrayUtil.transpose(ptsOut),
//            " ", ";" 
//				));
		
	}
	
	public static void test2d() {

		String fn = "/groups/jain/home/bogovicj/learning/advanced-imglib2/images/bee-1.tif";
		
		ImagePlus ip = IJ.openImage(fn);
		ImagePlusImg< ByteType, ? > img = ImagePlusImgs.from(ip);
		
		int[] dims = new int[]{	(int)img.dimension(0),
							   	(int)img.dimension(1)};
		
		System.out.println("img sz: " + ArrayUtil.printArray(dims));
		
//		int sPerDim = 9;
//		int N = numFacePts2d(sPerDim);
//		float[][] dest = new float[2][N];
//		
//		System.out.println("N: " + N);
//		
//		int k = 0;
//		facePts2dFloat(dims, sPerDim, dest, k );
//		System.out.println("face pts: \n" + ArrayUtil.printArray(dest));
		
		SimCrack sc = new SimCrack();
		sc.setDims(dims);
		
		int startX 		= 30;
		int crackLength = dims[0] - 30;
		int step        = 10;
		int fixedY 		= dims[1]/2;
		double crackHalfWidth = 10.0;
		
		sc.genSimpleCrackX2d( dims[0], startX, crackLength, step, fixedY, crackHalfWidth );
		sc.buildXfm2d();
		
		
		ImagePlus dest = new ImagePlus("dest", ip.getStack());
		System.out.println("in size:  " + ip.getWidth() + " " +  ip.getHeight() + " " + ip.getStackSize());
		System.out.println("dest size:  " + dest.getWidth() + " " +  dest.getHeight() + " " + dest.getStackSize());
		
		ThinPlateSplineMapping.mapInterval(sc.tps, ip, dest);
		
		String basename = "/groups/jain/home/bogovicj/projects/crackPatching/toyData/";
		
		IJ.save(dest, basename + "beeImg_crackX.tif");
		
		
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("starting\n");

//		test3d();
		test2d();
		
		
		
		System.out.println("\nfinished");
	}



}
