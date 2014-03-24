package net.imglib2.realtransform.spline;

import net.imglib2.RealLocalizable;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;


import java.util.ArrayList;

/**
 * Subclass of {@link KernelTransform} that implements a Thin-plate spline transformation.
 * Ported from itk's itkThinPlateSplineKernelTransform.hxx
 * <p>
 * M. H. Davis, a Khotanzad, D. P. Flamig, and S. E. Harms, 
 * “A physics-based coordinate transformation for 3-D image matching.,” 
 * IEEE Trans. Med. Imaging, vol. 16, no. 3, pp. 317–28, Jun. 1997. 
 *
 * @author Kitware (itk)
 * @author John Bogovic
 *
 */
public class ThinPlateSplineKernelTransform extends KernelTransform {

   protected static Logger logger = LogManager.getLogger(ThinPlateSplineKernelTransform.class.getName());

   public ThinPlateSplineKernelTransform(){
      super();
   }
   
   public ThinPlateSplineKernelTransform(Img<FloatType> img,
         ArrayList<RealLocalizable> srcPtList,
         ArrayList<RealLocalizable> tgtPtList) {
      super(img, srcPtList, tgtPtList);
   }

   @Override
   public void computeG(double[] pt, DenseMatrix64F mtx) {

      double nrm = Math.sqrt(normSqrd(pt));

      gMatrix = new DenseMatrix64F(ndims, ndims);
      CommonOps.setIdentity(gMatrix);
      CommonOps.scale(nrm,gMatrix);

   }

   @Override
   public double[] computeDeformationContribution(double[] thispt) {

		double[] l1 = new double[ndims];
      double[] res = new double[ndims];

      for (int lnd = 0; lnd < nLandmarks; lnd++) {

         sourceLandmarks.get(lnd).localize(l1);
         double[] diff = subtract(thispt, l1);
         double nrm = Math.sqrt(normSqrd(diff));

			logger.debug("dMatrix size: " + dMatrix.getNumRows() + " x " + dMatrix.getNumCols());
         for (int d = 0; d < ndims; d++) {
            res[d] += nrm * dMatrix.get(d, lnd);
         }
      }
      return res;
   }


}
