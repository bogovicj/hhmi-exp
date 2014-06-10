package net.imglib2.algorithms.crack.exps;

import ij.IJ;

import java.io.File;

import net.imglib2.algorithms.crack.*;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;
import net.imglib2.view.Views;

public class EdgelClusteringTests {

	public static void testClusteringToy(){
//		String destDir = "/Users/bogovicj/Documents/projects/crackStitching/edgeClustering";
		String destDir = "/groups/saalfeld/home/bogovicj/projects/crackStitching/edgeClustering/toy";
		
		// simulate a volume
		int[] sz = new int[]{200,200};
		Img<FloatType> img = ImgOps.createGradientImgY( sz, new FloatType());
		
		// generate simulated data
		float[] startPt = new float[]{   0, 100 };	
		float[] endPt   = new float[]{ 200, 100 };
		int crackLength = 11;
		SimCrack sc = new SimCrack( sz );
		sc.genCrack2d( startPt, endPt, crackLength, 5f );
		
		Img<FloatType> crackMask = sc.getCrackMask(img.factory(), img.firstElement());
		Img<FloatType> crackMaskSm = img.factory().create(img, img.firstElement());
		try {
//			Gauss3.gauss( 2, Views.extendMirrorDouble(crackMask), crackMaskSm);
			Gauss3.gauss( 2, crackMask, crackMaskSm );
		} catch (IncompatibleTypeException e) {
			e.printStackTrace();
		}
		
		//
		ImgOps.writeFloat( crackMask, destDir + File.separator + "crackMask.tif");
		ImgOps.writeFloat( crackMaskSm, destDir + File.separator + "crackMaskSmooth.tif");
		
		// cluster edgels
		int[] patchSize = new int[]{ 9,9 };
		CrackCorrection<FloatType> cc = new CrackCorrection<FloatType>( img, crackMaskSm, patchSize);
		cc.computeEdgels();
		
		System.out.println(" num edgels - " + cc.getEdgels().size());
		
		EdgelMatching<FloatType> em = new EdgelMatching<FloatType>();
		em.setSearchType( EdgelMatching.SearchTypes.COUNT );
		em.setEdgelSearchCount(27);
		em.setEdgels( cc.getEdgels() );

//		EdgelClustering<FloatType> ec = new EdgelClustering<FloatType>( em );
		EdgelClusteringRansac<FloatType> ec = new EdgelClusteringRansac<FloatType>( em );
		ec.cluster();
		
		Img<FloatType> clusterImg = img.factory().create( img, img.firstElement() );
		ec.makeEdgelClusterImgMem(clusterImg);
		ImgOps.writeFloat( clusterImg, destDir + File.separator + "edgelClusterImg.tif");
		
	}
	
	public static void testClusteringReal(){
//		String destDir = "/Users/bogovicj/Documents/projects/crackStitching/edgeClustering";
		String destDir = "/groups/saalfeld/home/bogovicj/projects/crackStitching/edgeClustering/closeup";
		
		// load a volume
		int downSampleFactor = 4;
		String imgfn = "/data-ssd1/john/projects/crackSegmentation/groundTruth/closeup/img_ds"+downSampleFactor+".tif";
		String maskfn = "/data-ssd1/john/projects/crackSegmentation/groundTruth/closeup/labels_interp_smooth_ds"+downSampleFactor+".tif";
		Img<FloatType> img =  ImagePlusAdapter.convertFloat( IJ.openImage(imgfn) );
		Img<FloatType> mask = ImagePlusAdapter.convertFloat( IJ.openImage(maskfn) );
		
		// cluster edgels
		int[] patchSize = new int[]{ 9,9 };
		CrackCorrection<FloatType> cc = new CrackCorrection<FloatType>( img, mask, patchSize);
		cc.computeEdgels();
		
		System.out.println(" num edgels - " + cc.getEdgels().size());
		
		EdgelMatching<FloatType> em = new EdgelMatching<FloatType>();
		em.setSearchType( EdgelMatching.SearchTypes.COUNT );
		em.setEdgelSearchCount(27);
		em.setEdgels( cc.getEdgels() );

//		EdgelClustering<FloatType> ec = new EdgelClustering<FloatType>( em );
		EdgelClusteringRansac<FloatType> ec = new EdgelClusteringRansac<FloatType>( em );
		ec.cluster();
		
		Img<FloatType> clusterImg = img.factory().create( img, img.firstElement() );
		ec.makeEdgelClusterImgMem(clusterImg);
		ImgOps.writeFloat( clusterImg, destDir + File.separator + "edgelClusterImg.tif");
		
	}
	
	public static void main(String[] args) {
		System.out.println("start");
		
//		testClusteringToy();
		testClusteringReal();
		
		System.out.println("finished");
		System.exit(0);
	}

}
