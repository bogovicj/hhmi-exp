package net.imglib2.algorithms.mgdm;

import java.util.*;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import net.imglib2.Cursor;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.integer.AbstractIntegerType;
import net.imglib2.type.numeric.real.AbstractRealType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.type.numeric.real.FloatType;
import edu.jhmi.rad.medic.libraries.ObjectProcessing;
import edu.jhmi.rad.medic.structures.BinaryHeap2D;
import edu.jhmi.rad.medic.utilities.Numerics;
import edu.jhu.ece.iacl.utility.ArrayReshape;
import edu.jhu.ece.iacl.utility.ArrayUtil;

/**
 * This class implements the MGDM level set decomposition for 2D and 3D images.
 * 
 * Given a discrete labeling, a "joint" fast marching simultaneously computes the
 * nearest neighbors and distance functions for each object.
 * 
 * @author John Bogovic
 * @version $Revision: $
 * 
 */
public class MgdmDecomposition <L extends AbstractIntegerType<L>> {
	
	private final static Logger logger = LogManager.getLogger(MgdmDecomposition.class.getName());
	
	public final static float UNKNOWN = -1f;
	public final static int   EMPTY = -1;
	
	protected final static double SQR2 = Math.sqrt(2.0);
	protected final static double SQR3 = Math.sqrt(3.0);
	
	protected final static double INF = 1e9;
	
	protected int[]   	dimSize;		// size of each spatial dimension
	protected float[] 	dimRes;			// grid spacing for each spatial dimension
	protected int		ndims;			// number of spatial dimensions
	protected int 		N; 				// total number of spatial points
	protected Img<L>    mask; 			// masking regions not used in computations
//	protected float[] 	res;			// repeated for fast marching
	
	protected Img<L> segmentation; 		// MGDM's segmentation
	protected Img<L> mgdmlabels; 		// MGDM's label maps
	protected int[] otherlabels; 		// MGDM's segmentation
	protected Img<FloatType> mgdmfunctions; // MGDM's pseudo level set mgdmfunctions
	protected float[]   tempObjLvlset; // temporary single-object levelset
	
	protected Img<DoubleType> phi;
	
	protected int nmgdm=3;	// number of mgdm functions to store
	protected int nobj;   // number of objects
	
	public double[]  res; 	// image resolutions
	

	protected float[] resoff;
	
	protected int[] objLabel; 	// list of original labels
	
	protected BinaryHeap2D heap; // the heap used in fast marching
	
	protected boolean areOriginalLabels;   	// keeps track of whether the current arrays
											// use the original labels or are mapped to indices
	
	
	// for debug 
	protected int 		debugX = -1;
	protected int 		debugY = -1;
	protected int 		debugZ = -1;
	protected int 		debugXYZ;
	protected boolean 	debug = true;
	

	
	
	/**
	 * Constructor for a Mgdm decomposition
	 * 
	 * @param segvol initial segmentation
	 * @param nmgdm the number of MGDM distance functions to store
	 * @param maskIn a mask outside of which MGDM computations will not be performed
	 */	
	public MgdmDecomposition(Img<L> segvol, int nmgdm, Img<L> maskvol){	
		
		segmentation = segvol.copy();
		
		// reshape and set the mask
		if(maskvol!=null){
			mask = maskvol.copy();
		}

		this.nmgdm=nmgdm;
		
		fastMarchingInitializationFromSegmentation( );
	}	

	/**
	 * Constructor for a Mgdm decomposition
	 * 
	 * @param segvol initial segmentation
	 * @param nmgdm the number of MGDM distance functions to store
	 */
	public MgdmDecomposition(Img<L> segvol, int nmgdm){
		this(segvol,nmgdm,null);
	}
	
	/**
	 * Constructor for a Mgdm decomposition/evolution.
	 * Defaults to using as many distance functions as 
	 * dimensions of the input image.
	 * 
	 * @param seg initial segmentation
	 */
	public MgdmDecomposition(Img<L> segvol){
		this(segvol, segvol.numDimensions());
	}
	
	private void initFunctions(Img<L> segvol){
		ndims = segvol.numDimensions();
		
		long[] dim = new long[ndims + 1];
		
		for(int d=0; d<ndims; d++){
			dim[d] = segvol.dimension(d);
		}
		dim[ndims] = nmgdm;
				 
		if(segmentation==null){
			segmentation = segvol.copy();
		}
		
		if(mgdmlabels==null){
			mgdmlabels = segvol.factory().create(dim, segvol.firstElement());
		}
		
		if(mgdmfunctions==null){
			ArrayImgFactory<FloatType> factory = new ArrayImgFactory<FloatType>();
			mgdmfunctions = factory.create(mgdmlabels, new FloatType());
		}
		
		if(phi==null){
			int[] sz = new int[ndims];
			ArrayUtil.fill(sz, 3);
			
			ArrayImgFactory<DoubleType> dfactory = new ArrayImgFactory<DoubleType>();
			phi = dfactory.create(sz, new DoubleType());
		}
		
		if(otherlabels==null){
			otherlabels = new int[N];
		}
	}
	
	/**
	 * Create a point at a definite location in a space of the dimensionality of
	 * the position.
	 *
	 * @param position
	 *            the initial position. The length of the array determines the
	 *            dimensionality of the space.
	 */
	public void setResolutions( final double... res )
	{
		if(res.length != ndims){
			logger.warn("length of arguments to setResolutions must equal number of dimensions = " + ndims);
			return;
		}
		
		for(int d=0; d<ndims; d++){
			this.res[d] = res[d];
		}
	}
	
	
	public final void toOriginalLabels() {
		
		if(areOriginalLabels){ // there's no work to do
			return;
		}
		
		Cursor<L> curs = segmentation.cursor();
		RandomAccess<L> mgdmRa = mgdmlabels.randomAccess();
		while ( curs.hasNext() ) {
			
			int i = curs.get().getInteger();
			
			if ( i < 0)
				continue;
			
			if( i >= objLabel.length){
				System.out.println("out of bounds - seg: " + i +"  objLab len: " + objLabel.length);
			}
			
			curs.get().setInteger(objLabel[i]);

			mgdmRa.setPosition(curs);
			mgdmRa.setPosition(0, ndims);
			for (int n = 0; n < nmgdm; n++) {
				int j = mgdmRa.get().getInteger();
				if ( j < 0){ continue; }
				mgdmRa.get().setInteger(  objLabel[j] );
			}
		}
		areOriginalLabels=true;
	}
	
	public int numDims(){
		return ndims;
	}
	
	public int getNumLabels(){
		return nobj;
	}
	
	/**
	 * Returns list of original labels
	 * @return
	 */
	public int[] getLabelList(){
		return objLabel;
	}
	
	
//	public float getDistance(int xyz, int n){
//		if(n<nmgdm){
//			return mgdmfunctions[n][xyz];
//		}else{
//			return Float.NaN;
//		}
//	}
	
//	public int getLabel(int xyz, int n){
//		if(n<nmgdm){
//			return mgdmlabels[n][xyz];
//		}else{
//			return -2;
//		}
//	}
//	
//	public int getLastNeighbor(int xyz){
//		return otherlabels[xyz];
//	}
//	
//	public void setSegmentation(int xyz, int value){
//		segmentation[xyz]=value;
//	}


	/**
	 * 
	 * Returns the signed distance function (phi) for an object, reconstructed
	 * from the MGDM decomposition
	 * 
	 * 
	 * @param pt the center point
	 * @param lb the object
	 * @return
	 */
	public final double labelSignedDistanceFunction(Localizable pt, int lb){
		
		RandomAccess<L> lbra = mgdmlabels.randomAccess();
		lbra.setPosition(pt);
		lbra.setPosition(0, ndims); //TODO double check me
		
		RandomAccess<FloatType> dra = mgdmfunctions.randomAccess();
		dra.setPosition(pt);
		dra.setPosition(0, ndims); //TODO double check me
		
		int l = lbra.get().getInteger();
		double out = 0;
		if( l == lb ){
			out = - dra.get().get();
		}
		
		for (int n = 0; n < nmgdm ; n++) {
			
			lbra.setPosition(n, ndims);
			l = lbra.get().getInteger();
			
			if ( l == lb){ break; }
			dra.setPosition(n, ndims);
			out += dra.get().get(); 
		}
		
		return out;
	}
	
	/**
	 * 
	 * Builds the signed distance function (phi) for an object, reconstructed
	 * from the MGDM decomposition in a neighborhood around the input {@link Localizable}.
	 * 
	 * 
	 * @param pt the center point
	 * @param lb the object
	 * @return
	 */
	public final void signedDistanceFunctionNeighborhood(Localizable pt, int lb){
		
		int[] pos = new int[ndims];
		pt.localize(pos);
		
		Cursor<DoubleType> c = phi.cursor();
		while(c.hasNext()){
			c.fwd();
			c.get().set( labelSignedDistanceFunction(c, lb));
		}
		
	}
	
	
//	public final int fastMarchingInitializationFromSegmentation(int[] init, float rx, float ry, float rz) {
//		return fastMarchingInitializationFromSegmentation(init,nmgdm,rx,ry,rz);
//	}
//	
//	public final int fastMarchingInitializationFromSegmentation(int[] init) {
//		return fastMarchingInitializationFromSegmentation(init,nmgdm,1f,1f,1f);
//	}
//	
//	public final int fastMarchingInitializationFromSegmentation(int[] init, int nmgdm) {
//		return fastMarchingInitializationFromSegmentation(init,nmgdm,1f,1f,1f);
//	}
	
	
	protected void fastMarchHeap(byte[] processed, double maxMarchDist ){
		
		float[] nbdist = new float[xoff.length];
		boolean[] nbflag = new boolean[xoff.length];
		
		float maxdist = 0.0f;
		// grow the labels and functions
		while (heap.isNotEmpty() && maxdist <= maxMarchDist ) {
			// extract point with minimum distance
			float dist = heap.getFirst();
			int xyz = heap.getFirstId();
			int lb = heap.getFirstState();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz] >= nmgdm)
				continue;

			// if there is already a label for this object, this is done
			boolean done = false;
			for (int n = 0; n < processed[xyz]; n++)
				if (mgdmlabels[n][xyz] == lb)
					done = true;
			if (done)
				continue;

			// update the distance functions at the current level
			mgdmfunctions[processed[xyz]][xyz] = dist;
			mgdmlabels[processed[xyz]][xyz] = lb;
			processed[xyz]++; // update the current level

			// keep track of distance if stopping at the narrow band
			maxdist = dist;

			// find new neighbors
			for (int k = 0; k < xoff.length; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];

				if (mask[xyzn]) {
					// must be in outside the object or its processed neighborhood
					boolean isprocessed = false;
					if (segmentation[xyzn] == lb)
						isprocessed = true;
					else {
						for (int n = 0; n < processed[xyzn]; n++)
							if (mgdmlabels[n][xyzn] == lb)
								isprocessed = true;
					}

					if (!isprocessed) {
						// compute new distance based on processed neighbors for
						// the same object
						for (int l = 0; l < xoff.length; l++) {
							
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							
							int xyznb = xyzn + xoff[l] + yoff[l] + zoff[l];
							// note that there is at most one value used here
							for (int n = 0; n < processed[xyznb]; n++) {
								if (mask[xyznb])
									if (mgdmlabels[n][xyznb] == lb) {
										nbdist[l] = mgdmfunctions[n][xyznb];
										nbflag[l] = true;
									}
							}
						}
						

						float newdist = minimumMarchingDistance(nbdist, nbflag, resoff);

						// add to the heap
						heap.addValue(newdist, xyzn, lb);
					}
				}
			}
		}

		// to create the MGDM functions, we need to copy the segmentation,forget
		// the last labels
		// and compute differences between distance functions
		// if(debug) System.out.println("transform into MGDM functions\n");

		for (int xyz = 0; xyz < N; xyz++)
			if (mask[xyz]) {
				// distance function difference
				for (int n = nmgdm - 1; n > 0; n--) {
					if (mgdmlabels[n][xyz] != EMPTY) {
						mgdmfunctions[n][xyz] = mgdmfunctions[n][xyz] - mgdmfunctions[n - 1][xyz];
					}
				}

				// label permutation
				otherlabels[xyz] = mgdmlabels[nmgdm - 1][xyz];
				for (int n = nmgdm - 1; n > 0; n--) {
					mgdmlabels[n][xyz] = mgdmlabels[n - 1][xyz];
				}
				mgdmlabels[0][xyz] = segmentation[xyz];

			}

	}

	/**
	 * Performs the joint fast marching on a discrete labeling, reshaped to 1D.
	 * @param init - the discrete labels
	 * @return zero if the fast marching was successful, non-zero if an error occured
	 */
	public final int fastMarchingInitializationFromSegmentation( ) {
		logger.info(" Fast martching over " + N + " spatial locations");
		
		

		return 0;
	}
	
	/**
	 * Performs the joint fast marching on a discrete labeling, reshaped to 1D.
	 * @param init - the discrete labels
	 * @return zero if the fast marching was successful, non-zero if an error occured
	 */
	public final int fastMarchingInitializationFromSegmentation_OLD( ) {

		
		if(nzPad==0){
			if(nxPad*nyPad != N ){
				System.err.println("Supplied dimensions (" +nx +" " + ny  + ") are incompatible with length of input label: " + N);
				return 1;
			}
		}else{
			if(nxPad*nyPad*nzPad != N ){
				System.err.println("Supplied dimensions (" +nx +" " + ny + " " + nz  + ") are incompatible with length of input label: " + N);
				return 1;
			}
		}
		
		System.out.println("Number of spatial locations: " + N);
		
		if(nzPad==0){
			// find labels in image
			objLabel = ObjectProcessing.listOrderedLabels(init, nxPad, nyPad);
			nobj = objLabel.length;
		}else{
			objLabel = ObjectProcessing.listOrderedLabels(init, nxPad, nyPad, nzPad);
			nobj = objLabel.length;
		}
		
		System.out.println("Found " + nobj + " objects");
		System.out.println(Arrays.toString(objLabel));
		
		// initialize heap
		if(heap==null){
			heap=new BinaryHeap2D(N,-1);
		}

		// initialize mgdm functions
		initFunctions();
		initOffset();
		
		// initialize the quantities
		// this replaces the labels with their rank indices (useful for making
		// things fast)
		for (int xyz = 0; xyz < N; xyz++) {
			// mgdm functions
			for (int n = 0; n < nmgdm; n++) {
				mgdmfunctions[n][xyz] = UNKNOWN;
				mgdmlabels[n][xyz] = EMPTY;
			}
			// segmentation
			
			int nlb = EMPTY;

			nlb = Arrays.binarySearch(objLabel, init[xyz]);
//			for (int n = 0; n < nobj; n++) {
//				if (objLabel[n] == init[xyz]) {
//					nlb = n;
//					continue;
//				}
//			}

			segmentation[xyz] = nlb;
			otherlabels[xyz] = EMPTY;
		}
		areOriginalLabels = false;
		
		//int[] afterSegLabels = ObjectProcessing.listOrderedLabels(segmentation, nxPad, nyPad, nzPad );
		//System.out.println(Arrays.toString(afterSegLabels));

		// computation variables
		byte[] processed = new byte[N]; // note: using a byte

		
		// compute the neighboring labels and corresponding distance functions
		// (! not the MGDM functions !)
		
		heap.reset();
		// initialize the heap from boundaries
		for (int xyz = 0; xyz < N; xyz++){
			if (mask[xyz]) {
				
				processed[xyz] = 0;
				// search for boundaries
				for (int k = 0; k < xoff.length; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (segmentation[xyzn] != segmentation[xyz])
						if (mask[xyzn]) {
							// add to the heap
							heap.addValue(0.5f * resoff[k], xyzn,
									segmentation[xyz]);
						}
				}
			}
		}

		fastMarchHeap(processed, Double.POSITIVE_INFINITY);

		return 0;
	}
	
	
	public final void fastMarchingReinitialization(){
		
		resetIsosurfaceBoundary();
		
		// computation variables
		byte[] processed = new byte[N]; // note: using a byte

		// compute the neighboring labels and corresponding distance functions
		// (! not the MGDM functions !)
		heap.reset();

		// initialize the heap from boundaries
		for (int xyz = 0; xyz < N; xyz++)
			if (mask[xyz]) {
				processed[xyz] = 0;
				segmentation[xyz] = mgdmlabels[0][xyz];
				// search for boundaries
				for (int k = 0; k < xoff.length; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (mgdmlabels[0][xyzn] != mgdmlabels[0][xyz])
						if (mask[xyzn]) {

							// add to the heap with previous value
							heap.addValue(mgdmfunctions[0][xyzn], xyzn,
									mgdmlabels[0][xyz]);
						}
				}
			}

		for (int xyz = 0; xyz < N; xyz++) {
			if (mask[xyz]) {
				for (int n = 0; n < nmgdm; n++) {
					mgdmlabels[n][xyz] = EMPTY;
				}
				otherlabels[xyz] = EMPTY;
			}
		}

		fastMarchHeap(processed, Double.POSITIVE_INFINITY);

		return;
	}
	
	protected final void resetIsosurfaceBoundary() {

		float[] nbdist = new float[xoff.length];
		boolean[] nbflag = new boolean[xoff.length];
		boolean boundary;

		float[] tmp = new float[N];
		boolean[] processed = new boolean[N];
		for (int xyz = 0; xyz < N; xyz++)
			if (mask[xyz]) {

				boundary = false;
				for (int l = 0; l < xoff.length; l++) {
					nbdist[l] = UNKNOWN;
					nbflag[l] = false;

					int xyznb = xyz + xoff[l] + yoff[l] + zoff[l];
					if (mgdmlabels[0][xyznb] != mgdmlabels[0][xyz]
							&& mask[xyznb]) {
						// compute new distance based on processed neighbors for
						// the same object
						nbdist[l] = Numerics.abs(mgdmfunctions[0][xyznb]);
						nbflag[l] = true;
						boundary = true;
					}
				}
				if (boundary) {
					tmp[xyz] = isoSurfaceDistance(mgdmfunctions[0][xyz],nbdist, nbflag, resoff);
					processed[xyz] = true;
				}
			}
		// once all the new values are computed, copy into original GDM function
		// (sign is not important here)
		for (int xyz = 0; xyz <N; xyz++) {
			if (processed[xyz])
				mgdmfunctions[0][xyz] = tmp[xyz];
			else
				mgdmfunctions[0][xyz] = UNKNOWN;
		}

		return;
	}
	

	
	/**
	 * The iso-surface distance computation (assumes a 6D array with opposite
	 * coordinates stacked one after the other) (the input values are all
	 * positive, the flags are true only if the iso-surface crosses)
	 * 
	 * @param cur
	 * @param val
	 * @param flag
	 * @param res sample spacing in each dimension
	 */
	public static final float isoSurfaceDistance(float cur, float[] val,
			boolean[] flag, float[] res) {

		if (cur == 0)
			return 0;

		float s;
		double dist;
		float tmp;
		s = 0;
		dist = 0;

		for (int n = 0; n < val.length; n += 2) {
			float ressqr = res[n] * res[n];
			if (flag[n] && flag[n + 1]) { // if both are across the boundary
				// Take the largest distance (aka the boundary surface is closer
				// to CURRENT point)
				tmp = Numerics.max(val[n], val[n + 1]);
				s = cur / (cur + tmp);
				dist += 1.0 / (s * s * ressqr);
			} else if (flag[n]) {
				s = cur / (cur + val[n]); // Else, take the boundary point
				dist += 1.0 / (s * s * ressqr);
			} else if (flag[n + 1]) {
				s = cur / (cur + val[n + 1]);
				dist += 1.0 / (s * s * ressqr);
			}
		}
		// triangular (tetrahedral?) relationship of height in right triangles
		// gives correct distance
		tmp = (float) Math.sqrt(1.0 / dist);

		// The larger root
		return tmp;
	}
	
	/**
	 * the Fast marching distance computation (!assumes a 6D array with opposite
	 * coordinates stacked one after the other, and res[0]==res[1],
	 * res[2]==res[3], res[4]==res[5]) This method differs from the above in
	 * that it incorporates image resolution into the distance computation
	 * 
	 * Solves for the distance value at 'this' voxel, x that is consistent with
	 * the eikonal equation: \nabla T = 1;
	 * 
	 * this is done by solving for the larger root of: ((x-phi1)/resx)^2 +
	 * ((y-phi2)/resy)^2 + ((z-phi3)/resz)^2 = 1
	 * 
	 * @param val
	 * @param flag
	 * @param res sample spacing in each dimension
	 */
	public static final float minimumMarchingDistance(float[] val,
			boolean[] flag, float[] res) {

		float s, s2; // s = rx*a + ry*b + rz*c; s2 = rx*rx*a*a + ry*ry*b*b +
		// rz*rz*c*c
		float tmp;
		float count;
		s = 0;
		s2 = 0;
		count = 0;

		for (int n = 0; n < val.length; n += 2) {
			float ressqr = res[n] * res[n];
			if (flag[n] && flag[n + 1]) {
				tmp = Numerics.min(val[n], val[n + 1]); // Take the smaller one
				// if both are processed
				s += tmp / ressqr;
				s2 += tmp * tmp / ressqr;
				count += (1.0 / ressqr);
			} else if (flag[n]) {
				s += val[n] / ressqr; // Else, take the processed one
				s2 += val[n] * val[n] / ressqr;
				count += (1.0 / ressqr);
			} else if (flag[n + 1]) {
				s += val[n + 1] / ressqr;
				s2 += val[n + 1] * val[n + 1] / ressqr;
				count += (1.0 / ressqr);
			}
		}

		// count must be greater than zero since there must be at least one
		// processed pt in the neighbors
		tmp = (s + (float) Math.sqrt((double) (s * s - count * (s2 - 1.0f))))
				/ (count);

		// The larger root
		return tmp;
	}
	
//	public boolean isDebug(int xyz){
//		return (xyz==debugXYZ);
//	}
//	
//	public int getDebugXYZ(){
//		return debugXYZ;
//	}
//	
//	
//	public void setDebugPoint(int debugX, int debugY, int debugZ){
//		this.debugX=debugX+2;
//		this.debugY=debugY+2;
//		this.debugZ=debugZ+2;
//		
//		// init debug parameters
//		debugXYZ = ArrayReshape.coordsToIndex(this.debugX, this.debugY,this.debugZ, nxPad, nyPad, nzPad);
//		
//	}
	
//	public String printDecomposition(int xyz){
//		String out = "";
//		for(int n=0; n<nmgdm; n++){
//			out += ("Label/Distance ["+n+"]: " + mgdmlabels[n][xyz] + "\t" + mgdmfunctions[n][xyz] + "\n") ;
//		}
//		out+="Last neighbor: " + otherlabels[xyz];
//		
//		return out;
//	}
	
//	public MgdmDecompositionIterator iterator(){
//		return new MgdmDecompositionIterator();
//	}
	
	
	public void semiFinalize(){
		mask = null;
		if(heap!=null){
			heap.finalize();
			heap = null;
		}
		System.gc();
	}
	
	public void finalize(){
		mgdmfunctions=null;
		mgdmlabels=null;
		segmentation=null;
		if(heap!=null){
			heap.finalize();
			heap = null;
		}
	}

//	/**
//	 * An iterator for the MGDM decomposition.
//	 * 
//	 * @author John Bogovic
//	 *
//	 */
//	public class MgdmDecompositionIterator implements Iterator<Integer>{
//		
//		// for iterating
//		protected int iterIdx;
//		protected int nextIdx;
//		
//		@Override
//		public boolean hasNext() {
//			nextIdx = iterIdx + 1;
//
//			if(nextIdx==N){ return false; } 
//
//			while(!mask[nextIdx]){
//
//				nextIdx++;
//
//				if(nextIdx==N){ // return false - out of bounds
//					return false;
//				}
//			}
//
//			return true;
//		}
//
//		@Override
//		public Integer next() {
//
//			if(nextIdx>=0){
//				iterIdx = nextIdx;
//			}else if(hasNext()){
//				iterIdx = nextIdx;
//			}else{
//				return -1;
//			}
//			return iterIdx;
//		}
//
//		@Override
//		public void remove() {
//			// do nothing
//		}
//
//		public void iteratorReset(){
//			iterIdx=0;
//		}
//	}
}

