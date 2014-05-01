package net.imglib2.algorithms.moments;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.RealType;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.ejml.alg.dense.decomposition.eig.SymmetricQRAlgorithmDecomposition;
import org.ejml.data.Complex64F;
import org.ejml.data.DenseMatrix64F;

import edu.jhu.ece.iacl.utility.ArrayUtil;
import edu.jhu.ece.iacl.utility.SortWithIndices;

/**
 * Compute image moments, and orientations derived therefrom.
 * 
 * @author John Bogovic
 *
 */
public class ImgMoment {
	protected static Logger logger = LogManager.getLogger(ImgMoment.class.getName());

	protected int        nd;
	
	double m;
	double[] m1;	
	double[] m2;
	
	protected double[] 	evals;
	protected double[] 	evecs;

	public static <T extends RealType<T>> double centralMoment(Img<T> img, int moment, int d){

		double centroid = moment(img, 1, d);

		logger.debug("moment " + moment );
		logger.debug("dim:" + d );

		int nd = img.numDimensions();
		double mom = 0; 

		if( d >= nd){
			logger.warn("Dimension for moment must be less than image dimensionality ... returning 0");
			return mom;
		}

		double[] loc = new double[nd];
		Cursor<T> lc = img.localizingCursor();
		while(lc.hasNext()){
			lc.next();
			lc.localize(loc);	

			mom += Math.pow( (loc[d] - centroid) , (double)moment ) * lc.get().getRealDouble();

		}

		return mom;
	}


	public static <T extends RealType<T>> double moment0( IterableInterval<T> img ){
		logger.debug("moment 0 ");

		double mom = 0; 

		Cursor<T> lc = img.localizingCursor();
		while(lc.hasNext()){
			lc.next();
			mom += lc.get().getRealDouble();
		}

		return mom;
	}

	public static <T extends RealType<T>> double[] moment1( IterableInterval<T> img ){
		logger.debug("moments 1 all ");

		int nd = img.numDimensions();
		double[] mom = new double[nd]; 

		double[] loc = new double[nd];
		Cursor<T> lc = img.localizingCursor();
		while(lc.hasNext()){
			lc.next();
			lc.localize(loc);	
			for(int d=0; d<nd; d++){
				mom[d] +=  loc[d] * lc.get().getRealDouble();
			}
		}

		logger.debug("result: " + ArrayUtil.printArray(mom));

		return mom;
	}

	public static <T extends RealType<T>> double moment1( IterableInterval<T> img, int d ){
		logger.debug("moment 1 ");

		double mom = 0; 
		int nd = img.numDimensions();

		if( d >= nd){
			logger.warn("Dimension for moment must be less than image dimensionality ... returning 0");
			return mom;
		}
		double[] loc = new double[nd];
		Cursor<T> lc = img.localizingCursor();
		while(lc.hasNext()){
			lc.next();
			lc.localize(loc);	
			mom +=  loc[d] * lc.get().getRealDouble();
		}

		logger.debug("result: " + mom);

		return mom;
	}

	public static <T extends RealType<T>> double[] moment2( IterableInterval<T> img ){
		logger.debug("moments 2 all ");

		int nd = img.numDimensions();
		double[] mom = new double[nd*nd]; 

		double[] loc = new double[nd];
		Cursor<T> lc = img.localizingCursor();
		while(lc.hasNext()){
			lc.next();
			lc.localize(loc);	
			int k = 0;
			for(int i=0; i<nd; i++)for(int j=0; j<nd; j++){
				mom[k++] +=  loc[i] * loc[j] * lc.get().getRealDouble();
			}
		}

		logger.debug("result: " + ArrayUtil.printArray(mom));

		return mom;
	}

	public static <T extends RealType<T>> double moment2( IterableInterval<T> img, int d ){
		int nd = img.numDimensions();
		logger.debug("moment 2 central for d " + d + " for " + nd + " dims");

		double mom = 0; 

		if( d >= nd){
			logger.warn("Dimension for moment must be less than image dimensionality ... returning 0");
			return mom;
		}
		double[] loc = new double[nd];
		Cursor<T> lc = img.localizingCursor();
		while(lc.hasNext()){
			lc.next();
			lc.localize(loc);	
			mom +=  loc[d] * loc[d] * lc.get().getRealDouble();
		}

		logger.debug("result: " + mom);

		return mom;
	}

	public static <T extends RealType<T>> double centralMoment2( Img<T> img, int d ){
		double c = moment1(img,d);
		int nd = img.numDimensions();

		logger.debug("moment 2 central for " + nd + " dims");

		double mom = 0; 

		if( d >= nd){
			logger.warn("Dimension for moment must be less than image dimensionality ... returning 0");
			return mom;
		}
		double[] loc = new double[nd];
		Cursor<T> lc = img.localizingCursor();
		while(lc.hasNext()){
			lc.next();
			lc.localize(loc);	
			mom +=  (loc[d] - c) * (loc[d] - c) * lc.get().getRealDouble();
		}

		return mom;
	}

	public static <T extends RealType<T>> double moment(IterableInterval<T> img, int moment, int d){

		int nd = img.numDimensions();

		logger.debug("moment " + moment + " for " + nd + " dims" );
		logger.debug("dim:" + d );

		double mom = 0; 

		if( d >= nd){
			logger.warn("Dimension for moment must be less than image dimensionality ... returning 0");
			return mom;
		}

		double[] loc = new double[nd];
		Cursor<T> lc = img.localizingCursor();
		while(lc.hasNext()){
			lc.next();
			lc.localize(loc);	

			mom += Math.pow( loc[d] , (double)moment ) * lc.get().getRealDouble();

		}

		return mom;
	}

	public <T extends RealType<T>> double[] orientation(IterableInterval<T> img){
		nd = img.numDimensions();
		logger.debug("orientation for " + nd + " dims" );
		double[] cov = new double[nd*nd];

		m  = moment0(img);	
		double msqr = m*m;
		if(msqr < 1e-9){ 
			logger.warn("... returning 0");
			return cov; 
		}
		m1 = moment1(img);	
		m2 = moment2(img);	

		int k = 0;
		for(int i=0; i<nd; i++)	for(int j=0; j<nd; j++) {
			cov[k] = m2[k]/m - m1[i]*m1[j]/msqr;
			k++;
		}

		return cov;
	}


	/**
	 * Returns null unless orientationEvalsEvecs is called
	 * @return the orientation eigenvalues
	 */
	public double[] getEvals(){
		return evals;
	}
	/**
	 * Returns null unless orientationEvalsEvecs is called
	 * @return the orientation eigenvectors
	 */
	public double[] getEvecs(){
		return evecs;
	}
	/**
	 * Returns null unless orientation is called
	 * @return the zero-th moment
	 */
	public double getM0(){
		return m;
	}
	/**
	 * Returns null unless orientation is called
	 * @return the first moment
	 */
	public double[] getM1(){
		return m1;
	}
	/**
	 * Returns null unless orientation is called
	 * @return the second moment
	 */
	public double[] getM2(){
		return m2;
	}
	
	/** 
	 * Computes a 12-vector
	 * The first element is the largest eigenvalue
	 * The next 3 elements are its associated eigenvector
	 * 
	 * The next 4 elements are the second eval and its evec
	 * and similarly for the last 4 elements.
	 *
	 * @param or Orientation computed from 'orientation' function
	 * @return  
	 */
	public void orientationEvalsEvecs(double[] or){

		DenseMatrix64F U = new DenseMatrix64F(ArrayUtil.reshape2D(or,nd,nd,false));

		SymmetricQRAlgorithmDecomposition decompU = new SymmetricQRAlgorithmDecomposition(true);
		decompU.decompose(U);
		int numU = decompU.getNumberOfEigenvalues();
		double[] evals = new double[nd];
		for(int i=0; i<numU; i++){ evals[i] = magsqr(decompU.getEigenvalue(i)); }

		int[] evalUinds = SortWithIndices.getIndexArray(nd);
		SortWithIndices.quicksort(evals, evalUinds);
		
		DenseMatrix64F Up  = permuteEvects(decompU, evalUinds);
		evecs = Up.getData();
		logger.debug(" " + Up );
		

		/* Older code using differently shaped output */

		//for(int evalIdx=0; evalIdx<3; evalIdx++){
		//	// eigenvectors are columns of the Up matrix
		//	for(int n=0; n<3; n++){
		//		out[n][evalIdx] = Up.getIndex(n, evalIdx);
		//	}   
		//}   
		//return out;

	}   

	public static DenseMatrix64F permuteEvects(SymmetricQRAlgorithmDecomposition in, int[] permutation){
		int N = in.getNumberOfEigenvalues();
		if(N==0){ return null; }

		DenseMatrix64F out = new DenseMatrix64F(N,in.getEigenVector(0).numRows);
		for(int i=0; i<N; i++){
			DenseMatrix64F evec = in.getEigenVector(permutation[out.numCols - i - 1]);
			logger.debug("\tevec: " + evec);
			for(int j=0; j<out.numCols; j++){
				out.set(i,j,evec.get(j));
			}
		}
		return out;
	}

	public static DenseMatrix64F permuteEvects(SymmetricQRAlgorithmDecomposition in, int[] permutation,float[] evals){
		int N = in.getNumberOfEigenvalues();
		if(N==0){ return null; }

		DenseMatrix64F out = new DenseMatrix64F(N,in.getEigenVector(0).numRows);
		for(int i=0; i<N; i++){
			DenseMatrix64F evec = in.getEigenVector(permutation[out.numCols - i - 1]);
			logger.debug("\tevec: " + evec);
			for(int j=0; j<out.numCols; j++){
				out.set(i,j, evals[i]*evec.get(j) );
			}
		}
		return out;
	}
	   

	public static double magsqr(Complex64F v){
		return (v.real*v.real + v.imaginary*v.imaginary);
	}

}
