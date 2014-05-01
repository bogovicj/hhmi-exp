package net.imglib2.algorithms.registration;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealRandomAccessible;
import net.imglib2.interpolation.InterpolatorFactory;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.InverseRealTransform;
import net.imglib2.realtransform.RealTransformRandomAccessible;
import net.imglib2.realtransform.RealViews;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import edu.jhu.ece.iacl.utility.ArrayUtil;

public class TransformTools {
	
	Logger logger = LogManager.getLogger( TransformTools.class.getName() );

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

	public static <T extends RealType<T>> RealTransformRandomAccessible<T, InverseRealTransform> xfmToView(
			AffineTransform3D xfm, 
			RandomAccessibleInterval<T> src, 
			InterpolatorFactory<T,
			RandomAccessible<T>> interpFactory) 
	{
		RealRandomAccessible<T> interpolant = Views.interpolate(
				Views.extendMirrorSingle(src), interpFactory);

		RealTransformRandomAccessible<T, InverseRealTransform> rv = 
			RealViews.transform( interpolant, xfm.inverse() );

		return rv;
	}
	
}
