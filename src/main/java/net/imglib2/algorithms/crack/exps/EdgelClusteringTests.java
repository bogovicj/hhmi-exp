package net.imglib2.algorithms.crack.exps;

import java.io.File;

import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.algorithms.crack.SimCrack;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;

public class EdgelClusteringTests {

	public static void testClusteringToy(){
		String destDir = "/Users/bogovicj/Documents/projects/crackStitching/edgeClustering";
		
		// simulate a volume
		int[] sz = new int[]{200,200};
		Img<FloatType> img = ImgOps.createGradientImgY( sz, new FloatType());
		
		// generate simulated data
		float[] startPt = new float[]{   5, 100 };	
		float[] endPt   = new float[]{ 195, 100 };
		int crackLength = 50;
		SimCrack sc = new SimCrack( sz );
		sc.genCrack2d( startPt, endPt, crackLength, 5f );
		
		Img<FloatType> crackMask = sc.getCrackMask(img.factory(), img.firstElement());
		Img<FloatType> crackMaskSm = img.factory().create(img, img.firstElement());
		try {
			Gauss3.gauss( 2, crackMask, crackMaskSm);
		} catch (IncompatibleTypeException e) {
			e.printStackTrace();
		}
		
		//
		ImgOps.writeFloat( crackMask, destDir + File.separator + "crackMask.tif");
		ImgOps.writeFloat( crackMaskSm, destDir + File.separator + "crackMaskSmooth.tif");
		
		// cluster edgels
		
	}
	
	public static void main(String[] args) {
		System.out.println("start");
		
		testClusteringToy();
		
		System.out.println("finished");
		System.exit(0);
	}

}
