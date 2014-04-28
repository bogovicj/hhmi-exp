package net.imglib2.algorithms.edge;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.ejml.data.DenseMatrix64F;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import edu.jhu.ece.iacl.utility.ArrayUtil;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithm.gradient.PartialDerivative;
import net.imglib2.algorithm.region.localneighborhood.Neighborhood;
import net.imglib2.algorithms.patch.PatchTools;
import net.imglib2.algorithms.region.localneighborhood.CrossShape;
import net.imglib2.algorithms.region.localneighborhood.CrossShape.NeighborhoodsAccessible;
import net.imglib2.img.Img;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.type.numeric.RealType;
import net.imglib2.*;
import net.imglib2.realtransform.*;
import net.imglib2.view.*;
import net.imglib2.view.composite.*;

public class EdgelTools {


	protected static Logger logger = LogManager.getLogger(EdgelTools.class
			.getName());


	public static <T extends RealType<T>> void laplacian( 
			RandomAccessible<T> src,
			Img<T> dest ){
		
		long[] gradDims = new long[dest.numDimensions() + 1];
		for ( int d=0; d<dest.numDimensions(); d++){
			gradDims[d] = dest.dimension(d);
		}
		gradDims[dest.numDimensions()] = dest.numDimensions();
		
		logger.debug("ndims: " + src.numDimensions());
		logger.debug("ndims: " + dest.numDimensions());
		
		// first partial derivatives
		Img<T> grad1 = dest.factory().create(gradDims, dest.firstElement());
		for ( int d=0; d<dest.numDimensions(); d++)
		{
			PartialDerivative.gradientCentralDifference(
					src, Views.hyperSlice( grad1, dest.numDimensions(), d), d);
		}
		
		// second partial derivatives
		Img<T> grad2 = dest.factory().create(gradDims, dest.firstElement());
		for ( int d=0; d<dest.numDimensions(); d++)
		{
			PartialDerivative.gradientCentralDifference(
					Views.extendMirrorDouble( Views.hyperSlice( grad1, dest.numDimensions(), d )), 
					Views.hyperSlice( grad2, dest.numDimensions(), d), d);
		}
		
		CompositeIntervalView<T, ? extends GenericComposite<T>> grad2VecView = Views.collapse(grad2);
		
		Cursor<T> dest_c = dest.cursor();
		CompositeIntervalView<T, ? extends GenericComposite<T>>.CompositeRandomAccess grad2Ra = grad2VecView.randomAccess();
		
		while( dest_c.hasNext() ) 
		{
			dest_c.fwd();
			grad2Ra.setPosition(dest_c);
			
			for ( long i=0; i<grad2.dimension(dest.numDimensions()); i++)
			{
				dest_c.get().add( grad2Ra.get().get(i) );
			}
		}
		
	}
	

	/**
	 * Finds the first zero-crossing of a 1d function with linear interpolation
	 * 
	 * @param src
	 */
	public static <T extends RealType<T>> double zeroXing1d(RandomAccessibleInterval<T> src)
	{
		double x = Double.NaN;
		
		if( src.numDimensions() != 1){
			logger.error("zeroXing1d accepts 1d intervals only - returning NaN");
			return x;
		}
		
		RandomAccess<T> ra = src.randomAccess();
		ra.setPosition(0, 0);
		double lastVal = ra.get().getRealDouble();
		
		
		for ( long i=src.min(0)+1; i<src.max(0); i++){
			ra.setPosition(i, 0);
			double curVal = ra.get().getRealDouble();
			
			if(curVal * lastVal < 0){ 
				// a sign change
				x = i - 1 - ( lastVal / ( curVal - lastVal )); 
				break;
			}
			else{
				// no sign change
				lastVal = curVal;
			}
		}
		
		if( Double.isNaN(x) ){
			logger.warn("Could not detect edge - returning NaN");
		}
			
		return x;
	}
	
			
	/**
	 * Find discrete local minima using a cross neighborhood.
	 * 
	 * @param src
	 * @param dest
	 */
	public static <T extends RealType<T>> void localAbsoluteMinNN(
			RandomAccessible<T> src,
			RandomAccessibleInterval<T> dest,
			double tol)
	{
		
		RandomAccess<T> srcRa  = src.randomAccess();
		RandomAccess<T> destRa = dest.randomAccess();
		
		CrossShape nbrhood = new CrossShape( 2, true);
		NeighborhoodsAccessible<T> windows = nbrhood.neighborhoods(dest);
		
		Cursor<Neighborhood<T>> winCurs = windows.cursor();
		while(winCurs.hasNext())
		{
			winCurs.fwd();
			
			srcRa.setPosition(winCurs);
			double ctrVal = Math.abs( srcRa.get().getRealDouble() );
			
			// check if any neighbors have a larger value
			boolean isMin = true;
			Cursor<T> nbrC = winCurs.get().cursor();
			while ( nbrC.hasNext() )
			{
				nbrC.fwd();
				
				srcRa.setPosition(nbrC);
				if( Math.abs(srcRa.get().getRealDouble()) + tol < ctrVal )
				{
					isMin = false;
					break;
				}
			}
			
			if( isMin )
			{
				destRa.setPosition(winCurs);
				destRa.get().setOne();
			}
			
		}
	}
	
	/**
	 * Replace this with a cross product if too slow..
	 * @param in an M x N matrix whose rows span a subspace of R^N
	 * @return an orthogonal basis for the subspace of R^N not spanned by in
	 */
	public static DenseMatrix64F remainderSubspace(DenseMatrix64F in){
		SimpleSVD<SimpleMatrix> svd = new SimpleSVD<SimpleMatrix>(in,false);
		SimpleMatrix nullspace = svd.nullSpace();
		return nullspace.getMatrix();
	}

	public static AffineTransform3D edgelTransformation(Edgel edgel){

		float[] nrm = edgel.getGradient();

		DenseMatrix64F normMtx = new DenseMatrix64F( new double[][]{ ArrayUtil.toDouble(nrm) } );
		DenseMatrix64F remMtx = remainderSubspace(normMtx);
		logger.debug("rem Mtx: " + remMtx );

		double[][] Ra = new double[3][4];
		for(int i=0; i<3; i++){
			for(int j=0; j<2; j++){
				Ra[i][j] = remMtx.get( i, j );
			}
			Ra[i][2] = normMtx.get(i);
		}

		logger.info(" Ra: \n" + ArrayUtil.printArray(Ra) );

		AffineTransform3D R = new AffineTransform3D();
		R.set(Ra);

		return R;
	}

	public static AffineTransform3D edgelToXfm(Edgel edgel, int[] midPt)
	{
		AffineTransform3D rotXfm = edgelTransformation(edgel);

		/******************************/
		AffineTransform3D xfm = rotXfm.copy();

		int ndims_in = midPt.length;
		double[] target = new double[ndims_in];


		xfm.apply( ArrayUtil.toDouble(midPt), target );
		double[] diff = ArrayUtil.subtract(
				ArrayUtil.toDouble(edgel.getPosition()), target);

		for (int i = 0; i < ndims_in; i++) {
			xfm.set(diff[i], i, ndims_in);
		}
		/******************************/

		return xfm;	
	}

	/**
	 * Returns a transformed view of the input source image relative to an edgel.
	 * The {@link Edgel} position will map to the midpoint the output patch.
	 * The +z axis of the output view corresponds to the gradient direction of 
	 * the edgel.
	 * 
	 * @param edgel the edgel
	 * @param src the source image
	 * @param patchSize the patch size
	 * @return the transformed view into the source image
	 */
	public static <T extends RealType<T>> RealTransformRandomAccessible<T, InverseRealTransform> edgelToView(Edgel edgel, RandomAccessibleInterval<T> src, int[] patchSize) 
	{
		logger.info(" edgel pos : " + ArrayUtil.printArray(edgel.getPosition()));
		logger.info(" edgel grad: " + ArrayUtil.printArray(edgel.getGradient()));
		logger.info(" edgel mag : " + edgel.getMagnitude());

		int[] midPt = PatchTools.patchSizeToMidpt( patchSize ); 	
		AffineTransform3D xfm = edgelToXfm(edgel, midPt);


		NLinearInterpolatorFactory<T> interpFactory = new NLinearInterpolatorFactory<T>();
//		NearestNeighborInterpolatorFactory<T> interpFactory = new NearestNeighborInterpolatorFactory<T>();

		RealRandomAccessible<T> interpolant = Views.interpolate(
				Views.extendMirrorSingle(src), interpFactory);

		RealTransformRandomAccessible<T, InverseRealTransform> rv = 
			RealViews.transform( interpolant, xfm.inverse() );

		return rv;
	}
}
