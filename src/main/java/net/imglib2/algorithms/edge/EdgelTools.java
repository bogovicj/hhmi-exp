package net.imglib2.algorithms.edge;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import org.ejml.data.DenseMatrix64F;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import edu.jhu.ece.iacl.utility.ArrayUtil;

import net.imglib2.RealRandomAccessible;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.img.Img;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.InverseRealTransform;
import net.imglib2.realtransform.RealTransformRandomAccessible;
import net.imglib2.realtransform.RealViews;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

public class EdgelTools {


	protected static Logger logger = LogManager.getLogger(EdgelTools.class
			.getName());


	/**
	 * Returns the integer coordinate of the midpoint of a patch
	 * with the specified size.  Is only accurate for patches
	 * with an odd size in every dimension.
	 * 
	 * @param patchSize size of the patch
	 * @return the midpoint coordinate
	 */
	public static int[] patchSizeToMidpt(int[] patchSize){ // determine translation

		int[] midPt = ArrayUtil.clone(patchSize);
		ArrayUtil.addInPlace(midPt, -1);
		ArrayUtil.divide(midPt, 2);

		return midPt;
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
	public static <T extends RealType<T>> RealTransformRandomAccessible<T,InverseRealTransform> edgelToView(Edgel edgel, Img<T> src, int[] patchSize) 
	{
		logger.info(" edgel pos : " + ArrayUtil.printArray(edgel.getPosition()));
		logger.info(" edgel grad: " + ArrayUtil.printArray(edgel.getGradient()));
		logger.info(" edgel mag : " + edgel.getMagnitude());

		int[] midPt = patchSizeToMidpt( patchSize ); 	
		AffineTransform3D xfm = edgelToXfm(edgel, midPt);


		NLinearInterpolatorFactory<T> interpFactory = new NLinearInterpolatorFactory<T>();

		RealRandomAccessible<T> interpolant = Views.interpolate(
				Views.extendMirrorSingle(src), interpFactory);

		RealTransformRandomAccessible<T, InverseRealTransform> rv = 
			RealViews.transform( interpolant, xfm.inverse() );

		return rv;
	}
}
