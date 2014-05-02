package net.imglib2.algorithms.registration;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealRandomAccessible;
import net.imglib2.algorithms.moments.ImgMoment;
import net.imglib2.interpolation.InterpolatorFactory;
import net.imglib2.realtransform.AbstractAffineTransform;
import net.imglib2.realtransform.AffineTransform;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.InverseRealTransform;
import net.imglib2.realtransform.InvertibleRealTransform;
import net.imglib2.realtransform.RealTransform;
import net.imglib2.realtransform.RealTransformRandomAccessible;
import net.imglib2.realtransform.RealViews;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import edu.jhu.ece.iacl.utility.ArrayUtil;

public class TransformTools {
	
	static Logger logger = LogManager.getLogger( TransformTools.class.getName() );

	/**
	 * Computes a rotation about an input center point
	 * 
	 * @param axis
	 * @param angle in radians
	 * @param center
	 * @return
	 */
	public static AffineTransform3D rotationCentered( int axis, double angle, double[] center)
	{
		/******************************/
		AffineTransform3D xfm = new AffineTransform3D();
		xfm.rotate( axis, angle );

		int ndims = center.length;
		double[] target = new double[ ndims ];
		
		xfm.apply( center, target );
		
		double[] diff = ArrayUtil.subtract( center, target );

		for (int i = 0; i < ndims; i++) {
			xfm.set(diff[i], i, ndims);
		}
		/******************************/

		return xfm;	
	}
	
	/**
	 * 
	 * @param R
	 * @param center
	 * @return
	 */
	public static DenseMatrix64F rotationCenteredMtx( DenseMatrix64F R, double[] center)
	{
		int ndims = center.length;
		
		logger.debug("R " + R);
		
		DenseMatrix64F affineXfm = new DenseMatrix64F( ndims + 1, ndims + 1 ); 
		DenseMatrix64F res = new DenseMatrix64F( ndims, 1 );
		
		DenseMatrix64F ctrVec = new DenseMatrix64F( ndims, 1 );
		
		ctrVec.setData( center );
		CommonOps.changeSign(ctrVec);
		
		CommonOps.mult(R, ctrVec, res);
		
		CommonOps.insert(  R, affineXfm, 0, 0);
		CommonOps.insert(res, affineXfm, 0, ndims);
		affineXfm.set( ndims, ndims, 1 );
		
		CommonOps.changeSign(ctrVec);
		for (int i = 0; i < ndims; i++) {
			affineXfm.set(i, ndims, 
					affineXfm.get(i,ndims) + ctrVec.get(i));
		}
		
		return affineXfm;
	}
	
	public static AffineTransform rotationCentered( DenseMatrix64F R, double[] center)
	{
		int ndims = center.length;
		AffineTransform xfm = toAffineTransform(R);
		double[] res = new double[ center.length ];
		xfm.apply( center, res);
		
		double[] diff = ArrayUtil.subtract( center, res );
		for (int i = 0; i < ndims; i++) {
			xfm.set(diff[i], i, ndims);
		}
		return xfm;
	}
	
	public static AffineTransform rotationCentered( double[][] R, double[] center)
	{
		int ndims = center.length;
		AffineTransform xfm = toAffineTransform(R);
		double[] res = new double[ center.length ];
		xfm.apply( center, res);
		
		double[] diff = ArrayUtil.subtract( center, res );
		for (int i = 0; i < ndims; i++) {
			xfm.set(diff[i], i, ndims);
		}
		return xfm;
	}

	public static <T extends RealType<T>> RealTransformRandomAccessible<T, InverseRealTransform> xfmToView(
			InvertibleRealTransform xfm, 
			RandomAccessible<T> src, 
			InterpolatorFactory<T, RandomAccessible<T>> interpFactory) 
	{
		RealRandomAccessible<T> interpolant = Views.interpolate(
				src, interpFactory);

		RealTransformRandomAccessible<T, InverseRealTransform> rv = 
			RealViews.transform( interpolant, xfm.inverse() );

		return rv;
	}
	
	/**
	 * Computes a rotation transformation from the source to target.
	 * @param src
	 * @param tgt
	 * @return
	 */
	public static <T extends RealType<T>> AffineTransform rotationPca(
			IterableInterval<T> src, IterableInterval<T> tgt ){
		
		int ndims = src.numDimensions();
		
		ImgMoment srcMom = new ImgMoment();
		double[] srcOr = srcMom.orientation(src);
		srcMom.orientationEvalsEvecs(srcOr);
		
		ImgMoment tgtMom = new ImgMoment();
		double[] tgtOr = tgtMom.orientation(tgt);
		tgtMom.orientationEvalsEvecs(tgtOr);
		
		DenseMatrix64F srcR = new DenseMatrix64F( ndims, ndims );
		srcR.setData( srcMom.getEvecs() );
		DenseMatrix64F srcMtx = rotationCenteredMtx(srcR, srcMom.centroid());
		
		DenseMatrix64F tgtR = new DenseMatrix64F( ndims, ndims );
		tgtR.setData( tgtMom.getEvecs() );
		DenseMatrix64F tgtMtx = rotationCenteredMtx(tgtR, tgtMom.centroid());
		
		DenseMatrix64F res =  new DenseMatrix64F( ndims+1, ndims+1 );
//		CommonOps.invert(srcMtx);
//		CommonOps.mult( tgtMtx, srcMtx, res );
		
		CommonOps.invert(tgtMtx);
		CommonOps.mult( srcMtx, tgtMtx, res );
		
		return toAffineTransform(res);
	}
	
	/**
	 * Computes a rotation transformation from the source to target.
	 * @param src
	 * @param tgt
	 * @return
	 */
	public static <T extends RealType<T>> AffineTransform rotationPca(
			IterableInterval<T> src, IterableInterval<T> tgt, double[] center ){
		
		int ndims = src.numDimensions();
		
		ImgMoment srcMom = new ImgMoment();
		double[] srcOr = srcMom.orientation(src);
		srcMom.orientationEvalsEvecs(srcOr);
		
		ImgMoment tgtMom = new ImgMoment();
		double[] tgtOr = tgtMom.orientation(tgt);
		tgtMom.orientationEvalsEvecs(tgtOr);
		
		DenseMatrix64F srcR = new DenseMatrix64F( ndims, ndims );
		srcR.setData( srcMom.getEvecs() );
		DenseMatrix64F srcMtx = rotationCenteredMtx(srcR, center);
		
		DenseMatrix64F tgtR = new DenseMatrix64F( ndims, ndims );
		tgtR.setData( tgtMom.getEvecs() );
		DenseMatrix64F tgtMtx = rotationCenteredMtx(tgtR, center);
		
		DenseMatrix64F res =  new DenseMatrix64F( ndims+1, ndims+1 );
//		CommonOps.invert(srcMtx);
//		CommonOps.mult( tgtMtx, srcMtx, res );
		
		CommonOps.invert(tgtMtx);
		CommonOps.mult( srcMtx, tgtMtx, res );
		
		return toAffineTransform(res);
	}
	
	public static AffineTransform toAffineTransform( DenseMatrix64F mtx ){
		int ndims = mtx.numRows - 1;
		AffineTransform xfm = new AffineTransform(ndims);
		for( int i=0; i<ndims; i++)for( int j=0; j<(ndims+1); j++){
			xfm.set( mtx.get(i, j), i, j);
		}
		return xfm;
	}
	public static AffineTransform toAffineTransform( double[][] mtx ){
		int ndims = mtx.length;
		AffineTransform xfm = new AffineTransform(ndims);
		for( int i=0; i<ndims; i++)for( int j=0; j<ndims; j++){
			xfm.set( mtx[i][j], i, j);
		}
		return xfm;
	}
	
	public static String printAffineTransform( AffineTransform xfm ){
		String s = "";
		int rows = xfm.numTargetDimensions();
		int cols = xfm.numSourceDimensions() + 1;
		for( int i=0; i<rows; i++){
			for( int j=0; j<cols; j++){
				s += " " + xfm.get(i, j);
			}
			s += "\n";
		}
		return s;
	}
}
