package net.imglib2.algorithms.crack;

import java.util.Arrays;
import java.util.Random;

import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import mpicbg.ij.ThinPlateSplineMapping;
import mpicbg.models.Point;
import mpicbg.models.PointMatch;
import mpicbg.models.TransformMesh;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.exception.ImgLibException;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.img.imageplus.ImagePlusImg;
import net.imglib2.img.imageplus.ImagePlusImgs;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.ByteType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import jitk.spline.KernelTransform;
import jitk.spline.KernelTransformFloatSeparable;
import jitk.spline.ThinPlateR2LogRSplineKernelTransform;
import jitk.spline.ThinPlateR2LogRSplineKernelTransformFloatSep;
import edu.jhu.ece.iacl.utility.ArrayUtil;
import edu.jhu.ece.iacl.utility.SimpleLinAlg;

//import mpicbg.models.TransformMesh;

public class SimCrack {

	int[] dims;
	int ndims;
	
	int[] 	dimsOut;
	int 	maxCoordOffset;
	
	ThinPlateR2LogRSplineKernelTransformFloatSep tps;
	
	int nCrackPts;
	float[][] 	crackCtrlPts;
	float[][] 	crackCtrlOffsets;
	
	CrackTransformMesh crackMesh;
	
	private int edgePtNumSamples = 21;
	float eps 			 = 0.1f;
	
	protected static Logger logger = LogManager.getLogger(SimCrack.class.getName());
	
	public SimCrack(){}
	
	public SimCrack(int[] dims )
	{
		setDims(dims);
	}
	
	public SimCrack(int[] dims, float[][] linePoints, float[][] offsets )
	{
		setDims(dims);
		this.crackCtrlPts = linePoints;
		this.crackCtrlOffsets = offsets;
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
	
	public void buildXfmMesh()
	{
		crackMesh = new CrackTransformMesh( dims );
		crackMesh.fromCrackParam( crackCtrlPts, crackCtrlOffsets );
		
	}
	
	public int[] outputSize( ){
		maxCoordOffset = -1;

		for( int i=0; i<crackCtrlOffsets.length; i++){
			for( int j=0; j<crackCtrlOffsets[0].length; j++){
				if ( crackCtrlOffsets[i][j] > maxCoordOffset ){
					maxCoordOffset = (int) Math.ceil( crackCtrlOffsets[i][j] );
				}
			}
		}
		maxCoordOffset += 2; // add a little extra for good measure

		int ndims = dims.length;
		dimsOut = new int[ ndims ];
		for( int d=0; d<ndims; d++){
			dimsOut[d] = (int)dims[d] + 2*maxCoordOffset;
		}
		return dimsOut;
	}
	
	public <T extends RealType<T>> Img<T> mapMesh( Img<T> src ){
		
		if ( dimsOut == null ){
			dimsOut = outputSize();
		}
		
		Img<T> dest = src.factory().create( dimsOut, src.firstElement() );
		IterableInterval<T> view = 
				Views.flatIterable(
					Views.offset( dest, maxCoordOffset, maxCoordOffset )
				);
		
		CrackTransformMeshMapping<CrackTransformMesh> map 
			= new CrackTransformMeshMapping<CrackTransformMesh>( crackMesh );
		
		map.mapInterpolated( src, view );
		
		return dest;
	}
	public <T extends NumericType<T>> Img<T> getCrackMask( ImgFactory<T> factory, T t )
	{
		Img<T> img = factory.create( dims, t );
		CrackTransformMeshMapping.mapMask( crackMesh, img );
		return img;
	}
	
	public void buildXfm2dSpline()
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
		tps.setLandmarks( tpsCtrlPtsSrc, tpsCtrlPtsTgt );
//		tps.setLandmarks( tpsCtrlPtsTgt, tpsCtrlPtsSrc);
		logger.info(" nLandmarks: " + tps.getNumLandmarks());
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
	public void genSimpleCrackX2d( int nx, int startX, int crackLength, int step, int fixedY, float crackHalfWidth )
	{
		crackCtrlPts 	 = new float[2][crackLength];
		crackCtrlOffsets = new float[2][crackLength];
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
	
	public void genSemiStepCrack( float[] startPt, float[] endPt, int crackLength,
									float crackWidth, float stepWidth, float maxSkew )
	{
		int stepIdx = crackLength / 2;
		System.out.println( "stepIdx " + stepIdx);
		
		int ndims = startPt.length;
		
		crackCtrlPts 	 = new float[crackLength][ndims];
		crackCtrlOffsets = new float[crackLength][ndims];
		
		crackCtrlPts[ 0 ] = startPt;
		crackCtrlPts[ crackLength-1 ] = endPt;
		
		float[] lineVec = ArrayUtil.subtract( endPt, startPt );
		ArrayUtil.normalizeLengthInPlace(lineVec);
		float[] lineVecPerp = SimpleLinAlg.orth2d( lineVec );
		
		float del = 1 / ( (float)crackLength);
		
		for( int i=1; i < (crackLength-1); i++)
		{
			float p = ( i + 1 ) * del;
			
			for( int d=0; d<ndims; d++){
				
				crackCtrlPts[ i ][ d ] = ( ( 1 - p ) * startPt[d] + (p * endPt[d]) );
				
				if( i > stepIdx )
				{
					crackCtrlPts[ i ][ d ] += stepWidth * lineVecPerp[d];
				}
			}	
		}
		
		crackCtrlOffsets[ 0 ] = ArrayUtil.multiply( lineVecPerp, eps);
		crackCtrlOffsets[ crackLength-1 ] = crackCtrlPts[ 0 ];
		float skewAmt = 0f;
		
		for( int i=1; i < (crackLength-1); i++)
		{
			if( i >= stepIdx )
			{
				skewAmt = maxSkew * ( ( crackLength - (float)i - 1) / ( crackLength - stepIdx - 1 ));
			}else {
				skewAmt = maxSkew * ( (float)i / ( stepIdx - 1 ));
			}
			
			System.out.println(" i " + i + "    skewAmt: " + skewAmt );
			
			
		}
	}
	
	public void genCrack2d( float[] startPt, float[] endPt, int crackLength,
			float distShp, float distScale, float skewMn, float skewVar)
	{
		genCrack2d( startPt, endPt, crackLength, distShp, distScale, skewMn, skewVar,
					new RandomDataGenerator());
	}
	
	public void genCrack2d( float[] startPt, float[] endPt, int crackLength,
							float distShp, float distScale, float skewMn, float skewVar,
							RandomDataGenerator rand ){
		
		int ndims = startPt.length;
		
		crackCtrlPts 	 = new float[crackLength][ndims];
		crackCtrlOffsets = new float[crackLength][ndims];
		
		crackCtrlPts[ 0 ] = startPt;
		crackCtrlPts[ crackLength-1 ] = endPt;
		
		double[] ptSpacing = new double[ crackLength - 2 ];
		for( int i=0; i<ptSpacing.length; i++){
			ptSpacing[i] = rand.nextUniform( 0, 1 );
		}
		Arrays.sort( ptSpacing );
		logger.debug(" pt Spacing \n" + ArrayUtil.printArray( ptSpacing ));
		
		for( int i=0; i<ptSpacing.length; i++){
			for( int d=0; d<startPt.length; d++){
				crackCtrlPts[ i+1 ][d] = (float)((1 - ptSpacing[i]) * startPt[d] + ptSpacing[i] * endPt[d]);  
			}	
		}
		
		float[] lineVec = ArrayUtil.subtract( endPt, startPt );
		ArrayUtil.normalizeLengthInPlace(lineVec);
		float[] lineVecPerp = SimpleLinAlg.orth2d( lineVec );
		
		for( int i=0; i<crackLength; i++)
		{
			double skewAmt = rand.nextGaussian( skewMn, skewVar);
			double dist = rand.nextGamma( distShp, distScale );
			
			for( int d=0; d<startPt.length; d++)
			{
				crackCtrlOffsets[i][d] = (float)( lineVec[d] * skewAmt + lineVecPerp[d] * ( 1 - skewAmt ) );
				crackCtrlOffsets[i][d] *= dist;
			}
		}
		
		logger.debug(" crack Pts:\n" + ArrayUtil.printArray( crackCtrlPts ));
		logger.debug("\n crack Offsets:\n" + ArrayUtil.printArray( crackCtrlOffsets ));
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
	public void genSimpleCrackXZ3d( int nx, int nz, int startX, int crackLength, int fixedY, float crackHalfWidth )
	{
		crackCtrlPts 	 = new float[3][crackLength*nz];
		crackCtrlOffsets = new float[3][crackLength*nz];
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
		
		float crackHalfWidth = 2.0f;
		sc.genSimpleCrackXZ3d(dims[0], dims[2], startX, crackLength, fixedY, crackHalfWidth);
		sc.buildXfm3d();
		
		System.out.println("sc: " + sc);
		
		FloatType b = new FloatType();
		Img<FloatType> img = ImgOps.createGradientImgY( dims, b);
		
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
		
		
		
		ThinPlateSplineMapping.mapInterval(sc.tps, ip.getProcessor(), dest.getProcessor());
		
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

		String basename = "/Users/bogovicj/Documents/projects/crackSim/grad1/";
		
//		String fn = "/groups/jain/home/bogovicj/learning/advanced-imglib2/images/bee-1.tif";
		String fn = basename + "grad.tif";
		
//		ImagePlus ip = IJ.openImage(fn);
//		ImagePlusImg< ByteType, ? > img = ImagePlusImgs.from(ip);
//		
//		int[] dims = new int[]{	(int)img.dimension(0),
//							   	(int)img.dimension(1)};
		
		int[] sz = new int[]{32,32}	;
		double[] w = new double[]{ 0.50, 0.50 };
		Img<FloatType> img = ImgOps.createGradientImg(sz, w, new FloatType());
		Img<FloatType> dest = img.factory().create( img, img.firstElement());
		
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
		sc.setDims(sz);
		
		int startX 		= 30;
		int crackLength = sz[0] - 15;
		int step        = 10;
		int fixedY 		= sz[1]/2;
		float crackHalfWidth = 3.0f;
		
		sc.genSimpleCrackX2d( sz[0], startX, crackLength, step, fixedY, crackHalfWidth );
		sc.buildXfm2dSpline();
		
		
//		ImagePlus dest = new ImagePlus("dest", ip.getStack());
//		System.out.println("in size:  " + ip.getWidth() + " " +  ip.getHeight() + " " + ip.getStackSize());
//		System.out.println("dest size:  " + dest.getWidth() + " " +  dest.getHeight() + " " + dest.getStackSize());
//		IJ.save(dest, basename + "beeImg_crackX.tif");
		
		mapInterval(sc.tps, img, dest);
		
		ImgOps.writeFloat(img, fn);
		ImgOps.writeFloat(dest, basename + "grad_crack.tif");
		
	}
	
	public static void test2dMesh(){
	
		String fnIn = "/Users/bogovicj/Documents/learning/advanced-imglib2/images/bee-1.tif";
		Img<FloatType> im = ImageJFunctions.convertFloat( IJ.openImage(fnIn));
		int[] sz = new int[]{ (int)im.dimension(0), (int)im.dimension(1) };
		String fnOut = "/Users/bogovicj/Documents/projects/crackSim/bee/bee_meshCrackRand_skew.tif";
		String fnMaskOut = "/Users/bogovicj/Documents/projects/crackSim/bee/bee_meshCrackRand_skew_mask.tif";
		
//		int[] sz = new int[]{ 200, 200 };
//		Img<FloatType> im = ImgOps.createGradientImgY( sz, new FloatType());
//		String fnOut = "/Users/bogovicj/Documents/projects/crackSim/grad1/grad_meshCrack_push.tif";
		
		SimCrack sc = new SimCrack( sz );

//		float[][] pts = new float[][]{ {100,100}, {115,100}, {135,100}, {150,100} };
//		float[][] off = new float[][]{ {0f,0.2f}, {0f,6f} ,   {0f,5f}, {0f,0.2f}  };
//		SimCrack sc = new SimCrack( sz, pts, off );
		
//		float[] startPt = new float[]{   0, 315 };
//		float[] endPt   = new float[]{ 600,   0 };
//		int crackLength =   25;
//		float distShp 	=   9f;
//		float distScale = 0.5f;
//		float skewMn 	=   0f;
//		float skewVar 	= 0.5f;
//		sc.genCrack2d(startPt, endPt, crackLength, distShp, distScale, skewMn, skewVar );
		
		float[] startPt = new float[]{   0, 315 };
		float[] endPt   = new float[]{ 600, 315 };
		int crackLength =   10;
		float crackWidth=   4f;
		float stepWidth =  12f;
		float maxSkew   =    5f;
		sc.genSemiStepCrack(startPt, endPt, crackLength, crackWidth, stepWidth, maxSkew);
		
//		System.out.println(" ctrlPts:\n " + ArrayUtil.printArray( sc.crackCtrlPts));
//		System.out.println(" ctrlOffsets:\n " + ArrayUtil.printArray( sc.crackCtrlOffsets));
		
//		sc.buildXfmMesh();
//		Img<FloatType> out = sc.mapMesh( im );
//		ImgOps.writeFloat( out, fnOut );
//		
//		Img<ByteType> mask = sc.getCrackMask( new ArrayImgFactory<ByteType>(), new ByteType());
//		ImgOps.writeByte( mask, fnMaskOut );
	}
	
	public final static <T extends NumericType<T>> void mapInterval(
			final KernelTransformFloatSeparable xfm,
			final Img<T> src, final Img<T> tgt )
	{
		NLinearInterpolatorFactory<T> interp = new NLinearInterpolatorFactory<T>();
		RealRandomAccess<T> sara = Views.interpolate( Views.extendZero(src), interp ).realRandomAccess();
		
		Cursor<T> tc = tgt.cursor();
		float[] pos = new float[src.numDimensions()]; 
		while( tc.hasNext() ){
			tc.fwd();
			tc.localize(pos);
			float[] srcPt  = xfm.transformPoint( pos );
			sara.setPosition( srcPt );
			tc.get().set( sara.get() );
			
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("starting\n");

//		test3d();
//		test2d();
		
		test2dMesh();
		
//		RandomDataGenerator r = new RandomDataGenerator();
//		for ( int i=0; i<50; i++){
//			System.out.println( r.nextUniform(0, 1));
//		}
		
		System.out.println("\nfinished");
	}



}
