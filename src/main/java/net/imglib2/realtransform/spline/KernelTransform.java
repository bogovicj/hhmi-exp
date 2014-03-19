package net.imglib2.realtransform.spline;

import java.util.ArrayList;
import java.util.Collection;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.ejml.data.*;
import org.ejml.factory.LinearSolver;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.ops.CommonOps;

import edu.jhu.ece.iacl.utility.ArrayUtil;
import net.imglib2.Dimensions;
import net.imglib2.RealLocalizable;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.type.numeric.real.FloatType;

/**
 * 
 * @author John Bogovic
 *
 */
public class KernelTransform {
	
	protected Dimensions dim;
   protected int ndims;

	protected DenseMatrix64F gMatrix;
	protected DenseMatrix64F pMatrix;
	protected DenseMatrix64F kMatrix;
	protected DenseMatrix64F dMatrix;
	protected DenseMatrix64F wMatrix;
	protected DenseMatrix64F lMatrix;
	protected DenseMatrix64F yMatrix;
	
	protected DenseMatrix64F I;

	protected double[][] aMatrix;
	protected double[] bVector;
	
	protected double 	stiffness = 0.25; // reasonable values take the range [0.0, 0.5]
	protected boolean	wMatrixComputeD = false; 
	
	protected ArrayList<RealLocalizable>  sourceLandmarks;
	protected ArrayList<RealLocalizable>  targetLandmarks;
	protected int 								nLandmarks;

	protected double[][] displacement;
	
	protected static Logger logger = LogManager.getLogger(KernelTransform.class.getName());
	
	public KernelTransform(){}
	
	public KernelTransform(Dimensions dim){
		logger.info("initializing");
		
		this.dim = dim;
		ndims = dim.numDimensions();

		sourceLandmarks = new ArrayList<RealLocalizable>();
		targetLandmarks = new ArrayList<RealLocalizable>();

		gMatrix = new DenseMatrix64F(ndims, ndims);

		I       = new DenseMatrix64F(ndims, ndims);
		for (int i=0; i<ndims; i++){
			I.set(i,i,1);
		}
		
	}
	public KernelTransform( Dimensions dim, Collection<RealLocalizable> srcPtsIn, Collection<RealLocalizable> tgtPtsIn ){
		this(dim);
		setLandmarks(srcPtsIn, tgtPtsIn);
	}

	public void setLandmarks(Collection<RealLocalizable> sourcePts, Collection<RealLocalizable> targetPts){
		if( sourcePts.size() != targetPts.size()){
			logger.error("Source and target landmark lists must be the same size");
			return;
		}
		setSourceLandmarks(sourcePts);
		setTargetLandmarks(targetPts);
		
		pMatrix = new DenseMatrix64F( (ndims * nLandmarks), ( ndims * (ndims + 1)) );
		dMatrix = new DenseMatrix64F( ndims, nLandmarks);
		wMatrix = new DenseMatrix64F( ndims, nLandmarks);
		kMatrix = new DenseMatrix64F( ndims * nLandmarks, ndims * nLandmarks);
		yMatrix = new DenseMatrix64F( ndims * ( nLandmarks + ndims + 1), 1 );
		lMatrix = new DenseMatrix64F( ndims * ( nLandmarks + ndims + 1),
									  ndims * ( nLandmarks + ndims + 1) );
		wMatrix = new DenseMatrix64F( (ndims * nLandmarks) + ndims * ( ndims + 1),
				  					  1 );
		
		aMatrix = new double[ndims][ndims];
		bVector = new double[ndims];
		displacement = new double[nLandmarks][ndims];
		
		
	}

	protected void printLandmarks(){
		double[] tmp = new double[ndims];

		System.out.println("\nSOURCE LANDMARKS:");
		for( RealLocalizable lm : sourceLandmarks){
			lm.localize(tmp);
			System.out.println("  " + ArrayUtil.printArray(tmp));
		}
	
		System.out.println("\nTARGET LANDMARKS:");
		for( RealLocalizable lm : targetLandmarks){
			lm.localize(tmp);
			System.out.println("  " + ArrayUtil.printArray(tmp));
		}
		System.out.println("\n");
	}

	protected void setSourceLandmarks(Collection<RealLocalizable> sourcePts){
		for (RealLocalizable l : sourcePts ){
			sourceLandmarks.add(l);
		}
		nLandmarks = sourceLandmarks.size();
	}
	
	protected void setTargetLandmarks(Collection<RealLocalizable> targetPts){
		for (RealLocalizable l : targetPts ){
			targetLandmarks.add(l);
		}
	}
	
	public double[] transformPoint(double[] pt){

      double[] result = new double[ndims];
      // affine part
		for(int i=0; i<ndims; i++){
			for(int j=0; j<ndims; j++){
            result[i] += aMatrix[i][j] * pt[j];
         }
      }

      // translational part
      for(int i=0; i<ndims; i++){
         result[i] += bVector[i] + pt[i];
      }

		return result;
	}
	
	public DenseMatrix64F computeReflexiveG(){
		CommonOps.fill(gMatrix, 0);
		for (int i=0; i<ndims; i++){
			gMatrix.set(i,i, stiffness);
		}
		return gMatrix;
	}

	public double[] computeDeformationContribution( double[] thispt ){

		double[] result = new double[ndims];
		double[] l1 = new double[ndims];

		for( int lnd=0; lnd<nLandmarks; lnd++){
			sourceLandmarks.get(lnd).localize(l1);
			computeG( result, gMatrix );
			for (int i=0; i<ndims; i++) for (int j=0; j<ndims; j++){
				result[j] += gMatrix.get(i,j) * dMatrix.get(i,lnd);
			}

		}
		return result;
	}

	public void computeD(){

		for( int i=0; i<nLandmarks; i++ ){
			displacement[i] = subtract(	targetLandmarks.get(i),
													sourceLandmarks.get(i));
		}
	}	

	/**
	 *
	 */
	protected double[] subtract(RealLocalizable p1, RealLocalizable p2){
		int nd = p1.numDimensions();
		double[] out = new double[nd];
		for (int d=0; d<nd; d++){
			out[d] = p1.getDoublePosition(d) - p2.getDoublePosition(d);
		}
		return out;
	}
 
	protected double normSqrd( double[] v ){
     double nrm = 0;
      for(int i=0; i<v.length; i++){
         nrm += v[i]*v[i]; 
      }
      return nrm;
   }

	protected double[] subtract(double[] p1, double[] p2){
		int nd = p1.length; 
		double[] out = new double[nd];
		for (int d=0; d<nd; d++){
			out[d] = p1[d] - p2[d];
		}
		return out;
	}
	
	protected double[] subtract(double[] p1, double[] p2, double[] out){
		int nd = out.length; 
		for (int d=0; d<nd; d++){
			out[d] = p1[d] - p2[d];
		}
		return out;
	}

	public void computeW(){
		computeL();
		computeY();

		// solve linear system 
		LinearSolver<DenseMatrix64F> solver =  LinearSolverFactory.pseudoInverse(true);
		solver.setA(lMatrix);
		solver.solve(yMatrix, wMatrix);

		reorganizeW();
	}

	public void computeL(){

		computeP();
		computeK();

		CommonOps.insert( kMatrix, lMatrix, 0, 0 );
		CommonOps.insert( pMatrix, lMatrix, 0, kMatrix.getNumCols() );
		CommonOps.transpose(pMatrix);

		CommonOps.insert( pMatrix, lMatrix, kMatrix.getNumRows(), 0 );
		CommonOps.insert( kMatrix, lMatrix, 0, 0 );
	
      // bottom left O2 is already zeros after initializing 'lMatrix'
      
	}

	public void computeP(){
		
		DenseMatrix64F tmp = new DenseMatrix64F(ndims,ndims);

		double[] p1  = new double[ndims];

		for( int i=0; i<nLandmarks; i++ ){

			sourceLandmarks.get(i).localize(p1);

			for( int j=0; j<ndims; j++ ){
				CommonOps.scale(p1[j], I, tmp);
				CommonOps.insert( tmp, pMatrix,  i*ndims, j*ndims );

			}
			CommonOps.insert( I, pMatrix,  i*ndims, ndims*ndims );
		}
	}

	public void computeK(){

	  computeD();

      double[] p1  = new double[ndims];
      double[] p2  = new double[ndims];
      double[] res = new double[ndims];

		for( int i=0; i<nLandmarks; i++ ){

         sourceLandmarks.get(i).localize(p1);

			DenseMatrix64F G = computeReflexiveG();

			CommonOps.insert(G, kMatrix, i * ndims, i * ndims);
			
			for( int j=i; j<nLandmarks; j++ ){

				sourceLandmarks.get(j).localize(p2);

				subtract( p1, p2, res );
				computeG(res, G);			

				CommonOps.insert(G, kMatrix, i * ndims, j * ndims);
				CommonOps.insert(G, kMatrix, j * ndims, i * ndims);

			}
		}
	}

	public void computeY(){

		CommonOps.fill( yMatrix, 0 );

		for (int i=0; i<nLandmarks; i++) {
			for (int j=0; j<ndims; j++) {
				yMatrix.set( i*ndims + j, 0, displacement[i][j]);
			}
		}
		for (int i=0; i< ndims * (ndims + 1); i++) {
			yMatrix.set( nLandmarks * ndims + i, 0, 0);
		}

	}


	public void reorganizeW(){

		
		int ci = 0;

		for( int lnd=0; lnd<nLandmarks; lnd++){
			for (int i=0; i<ndims; i++) {
				dMatrix.set(i, lnd, wMatrix.get(ci, 0));
				ci++;
			}	
		}

		for( int j=0; j<ndims; j++) for (int i=0; i<ndims; i++) {
			aMatrix[i][j] =  wMatrix.get(ci,0);
			ci++;
		}

		for( int k=0; k<ndims; k++) {
			bVector[k] = wMatrix.get(ci, 0);
			ci++;
		}

		wMatrix = null;
		
	}

	public double[] transformPoint(RealLocalizable loc){
		
		double[] pt 	 = new double[ndims];
		loc.localize(pt);

		double[] result = computeDeformationContribution( pt );

		for (int i = 0; i < ndims; i++)
			for (int j = 0; j < ndims; j++) {
				result[i] += aMatrix[i][j] * pt[j];
			}

		return result;
	}


	public void computeG( double[] pt, DenseMatrix64F mtx){  }

	public void Modified(){
		
	}

	public static void main(String[] args) {

		FloatType t = new FloatType(0f);
		long[] d = new long[]{20,20};
		ArrayImg<FloatType,FloatType> img = new ArrayImg<FloatType,FloatType>(t, d, 0);
		
		KernelTransform xfm = new KernelTransform(img); 
		System.out.println("xfm: " + xfm);
				

		// set landmarks
		
		// computeW

	}

}
