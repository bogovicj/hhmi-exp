package net.imglib2.education;

import ij.IJ;
import ij.ImagePlus;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import jitk.spline.ThinPlateR2LogRSplineKernelTransformFloat;
import mpicbg.ij.ThinPlateSplineMapping;
import mpicbg.trakem2.transform.ThinPlateSplineTransform;

public class TpsBeeExample {
	public static void testBoxDeformation()
	{
		

		int nx = 640;
		int ny = 374;
		
		int[] xr = new int[]{200, 400};
		int[] yr = new int[]{150, 230};
		
//		int[] dx = new int[]{-20, 0, 20};
//		int[] vx = new int[]{-40, 0, 40};
//		
//		int[] dy = new int[]{-20, 0, 20};
//		int[] vy = new int[]{-40, 0, 40};
		
		int[] dx = new int[]{ -20, - 5, 0 };
		int[] vx = new int[]{   0, -10, 0 };
		
		int[] dy = new int[]{ -20, - 5, 0 };
		int[] vy = new int[]{   0, -10, 0 };
		
//		int[] dy = new int[]{-20, 0, 20};
//		int[] vy = new int[]{-40, 0, 40};
				
		int npts = xr.length * yr.length * dx.length;
		
		
		float[][] srcPts = new float[2][npts+4];
		float[][] tgtPts = new float[2][npts+4];
		
		System.out.println("srcPts sz  " + srcPts.length + " x " + srcPts[0].length);
		
		int k = 0;
		for(int xi=0; xi<xr.length; xi++)for(int yi=0; yi<yr.length; yi++){

			// in here we're at one point

			for( int di=0; di<dx.length; di++)
			{
				
				srcPts[0][k] = xr[xi] + dx[di];
				srcPts[1][k] = yr[yi] + dy[di];
				
				tgtPts[0][k] = srcPts[0][k] + vx[di];
				tgtPts[1][k] = srcPts[1][k] + vy[di];
				
				k++;
			}
		}
		
		// add boundary clamps
		srcPts[0][k] = 0;
		tgtPts[0][k] = 0;
		srcPts[1][k] = 0;
		tgtPts[1][k] = 0;
		k++;
		
		srcPts[0][k] = nx;
		tgtPts[0][k] = nx;
		srcPts[1][k] = 0;
		tgtPts[1][k] = 0;
		k++;
		
		srcPts[0][k] = 0;
		tgtPts[0][k] = 0;
		srcPts[1][k] = ny;
		tgtPts[1][k] = ny;
		k++;
		
		srcPts[0][k] = nx;
		tgtPts[0][k] = nx;
		srcPts[1][k] = ny;
		tgtPts[1][k] = ny;
		k++;
		
		System.out.println(" solving for transformation ..." );
		
		String dstImgFn = "/Users/bogovicj/Documents/learning/advanced-imglib2/images/bee-1-tpsXfm2.tif";
		testTpsBee( srcPts, tgtPts, dstImgFn );

	}
	
	public static void testLegsDeformation(){
		int nx = 640;
		int ny = 374;
		
		float[][] srcPts = new float[][]{
				{0,  0, nx, nx, 171, 297, 444, 221, 293, 395},
				{0, ny, 0 , ny, 262, 265, 261, 211, 211, 211 }
		};
		
		float[][] tgtPts = new float[][]{
				{0,  0, nx, nx, 127, 297, 490, 221, 293, 395 },
				{0, ny,  0, ny, 313, 313, 313, 211, 211, 211 }
		};
		
		
		String dstImgFn = "/Users/bogovicj/Documents/learning/advanced-imglib2/images/bee-1-tpsXfmLEGS-REV.tif";
//		String dstXfmFn = "/Users/bogovicj/Documents/learning/advanced-imglib2/images/bee-1-tpsXfm.txt";
		testTpsBee( srcPts, tgtPts, dstImgFn );

	}
	
	public static void testTpsBee(float[][] srcPts, float[][] tgtPts, String dstImgFn ){
		
		String srcImgFn = "/Users/bogovicj/Documents/learning/advanced-imglib2/images/bee-1.tif";
		
		ImagePlus srcImg = IJ.openImage(srcImgFn);
		ImageProcessor srcIp = srcImg.getProcessor();
		int nx = srcIp.getWidth();
		int ny = srcIp.getHeight();
		
		ThinPlateR2LogRSplineKernelTransformFloat tps = new 
				ThinPlateR2LogRSplineKernelTransformFloat(2, tgtPts, srcPts);
//		ThinPlateR2LogRSplineKernelTransformFloat tps = new 
//				ThinPlateR2LogRSplineKernelTransformFloat(2, srcPts, tgtPts);
		
//		tps.setDoAffine(false);
		tps.computeW();
		
//		System.out.println( tps.getAffine() );
		
		ThinPlateSplineTransform xfm  = new ThinPlateSplineTransform(tps);
		ThinPlateSplineTransform xfm2 = new ThinPlateSplineTransform();
		
		String str = xfm.toDataString();
		System.out.println("\n" + str + "\n");
		
		xfm2.init(str);
		
		float[] loc = new float[]{127, 313};
		xfm.applyInPlace(loc);
		
		System.out.println(" loc: " + loc[0] + " " + loc[1]);
		
		ImageProcessor dstIp = new ByteProcessor( nx, ny );
		ThinPlateSplineMapping.mapInterval(tps, srcIp, dstIp);
		
		ImagePlus dstImg = new ImagePlus(srcImg.getTitle()+"_tps", dstIp);
		
		IJ.save(dstImg, dstImgFn);
		
	}
	
	
	public static void main(String[] args){
		System.out.println("Starting");
		
//		testBoxDeformation();
		testLegsDeformation();
		
		System.out.println("Finished");
		System.exit(0);
	}

}
