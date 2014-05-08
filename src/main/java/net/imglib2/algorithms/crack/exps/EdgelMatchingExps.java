package net.imglib2.algorithms.crack.exps;

import java.util.ArrayList;

import edu.jhu.ece.iacl.utility.ArrayUtil;
import ij.IJ;
import net.imglib2.RealPoint;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithm.edge.EdgelFloat;
import net.imglib2.algorithm.edge.SubpixelEdgelDetection;
import net.imglib2.algorithm.edge.SubpixelEdgelDetectionFloat;
import net.imglib2.algorithms.crack.CrackCorrection;
import net.imglib2.algorithms.crack.EdgelMatching;
import net.imglib2.collection.KDTree;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.neighborsearch.RadiusNeighborSearchOnKDTree;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;

public class EdgelMatchingExps {

	public static void tryMatchRegistration() 
	{
		String imgfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/crackVolDown_cp.tif";
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_smooth.tif";
		int[] patchSize = new int[] { 19, 19, 7 };
		
		Img<FloatType> img =  ImagePlusAdapter.convertFloat( IJ.openImage(imgfn) );
		Img<FloatType> mask = ImagePlusAdapter.convertFloat( IJ.openImage(maskfn) );

		CrackCorrection<FloatType> cc = new CrackCorrection<FloatType>(
				img, mask, patchSize);
		cc.computeEdgels();
		
		int i = cc.edgelIdxNearest(new double[]{67,290,13});
		int j = cc.edgelIdxNearest(new double[]{69,311,13});
		
		System.out.println(" i: " + i + "  " + cc.getEdgels().get(i));
		System.out.println(" j: " + j + "  " + cc.getEdgels().get(j));
	
		cc.debugOutDir = "/groups/jain/home/bogovicj/projects/crackPatching/edgelRegPatch";
//		cc.registerEdgelsOrient( cc.getEdgels().get(i), cc.getEdgels().get(j), i);
		
	}
	
	public static void tryMatching() 
	{
//		String imgfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/crackVolDown_cp.tif";
//		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_smooth.tif";
		
		int downSampleFactor = 4;
		double searchRadius = 100;
		
//		String imgfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/groundTruth/closeup/img.tif";
//		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/groundTruth/closeup/labels_interp_smooth.tif";
		
		String imgfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/groundTruth/closeup/img_ds"+downSampleFactor+".tif";
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/groundTruth/closeup/labels_interp_smooth_ds"+downSampleFactor+".tif";
		int[] patchSize = new int[] { 31, 31, 19 };
		
		double[] testPt = new double[]{ 159,253,73 }; 
		double[] resPt  = new double[]{ 170, 312, 136 };
		
		
		ArrayUtil.divide(testPt, downSampleFactor);
		ArrayUtil.divide(resPt, downSampleFactor);
		
		double[] diff   = ArrayUtil.clone(resPt);
		ArrayUtil.subtractInPlace(diff, testPt);
		double dist = Math.sqrt(ArrayUtil.sumSquares(diff));
		System.out.println( " dist " + dist);
		
		ArrayUtil.divide( testPt, (double)downSampleFactor);
		System.out.println( "testPt: " +  ArrayUtil.printArray(testPt));
		
		Img<FloatType> img =  ImagePlusAdapter.convertFloat( IJ.openImage(imgfn) );
		Img<FloatType> mask = ImagePlusAdapter.convertFloat( IJ.openImage(maskfn) );

		CrackCorrection<FloatType> cc = new CrackCorrection<FloatType>(
				img, mask, patchSize);
		
		cc.computeEdgels();
		
//		cc.edgelIdxImg();
		
//		ArrayList<Edgel> edgels = cc.getEdgels();
//		int i =0;
//		for ( Edgel e : edgels ){
//			System.out.println( ArrayUtil.printArray( e.getPosition() ));
//			
//			if(i>50) break;
//		}
		

		
		EdgelMatching<FloatType> em = new EdgelMatching<FloatType>(img, mask, patchSize);
		em.debugDir = "/groups/jain/home/bogovicj/projects/crackPatching/cropEdgelMatch";
		em.debugSuffix = "ds"+downSampleFactor;
//		em.debugSuffix = "";
		
//		em.setSearchType( EdgelMatching.SearchTypes.COUNT );
//		em.setEdgelSearchCount(3000);
		
		em.setEdgelSearchRadius( searchRadius / (downSampleFactor) );
		
		em.setEdgels(cc.getEdgels());
		
//		em.computeAllAffinities();
//		em.testAffinitiesReg();
		
		
		int testEdgelIdx = cc.edgelIdxNearest(testPt);
		em.debug_i  = testEdgelIdx;
		
		Edgel testEdgel = cc.getEdgels().get(testEdgelIdx);
		em.computeAffinities( testEdgel );
		
	}
	
	
	public static void stupdidRadiusSearchTest()
	{
		
		ArrayList<RealPoint> ptlist = new ArrayList<RealPoint>();
		
		for(double x=0; x<5; x++) for(double y=0; y<5; y++){
			ptlist.add(new RealPoint(x,y));
		}
		
		KDTree<RealPoint> tree = new KDTree<RealPoint>( ptlist, ptlist );
		RadiusNeighborSearchOnKDTree<RealPoint> search =
				new RadiusNeighborSearchOnKDTree<RealPoint>( tree );
		
		search.search( new RealPoint(0.5,0.5), 1, false);
		
		System.out.println( " found " + search.numNeighbors() + " matches");
		
	}
	
	public static void edgelTypeValidate() {
		
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_smooth.tif";

		
		ArrayImgFactory<FloatType> ffactory = new ArrayImgFactory<FloatType>();
		ArrayImgFactory<DoubleType> dfactory = new ArrayImgFactory<DoubleType>();
		
		Img<FloatType> mask = ImagePlusAdapter.convertFloat( IJ.openImage(maskfn) );
	
		float thresh = 40.0f;
		
		ArrayList<Edgel> 		edgels  = SubpixelEdgelDetection.getEdgels(mask, ffactory, thresh);
		ArrayList<EdgelFloat> 	edgelsf = SubpixelEdgelDetectionFloat.getEdgels(mask, ffactory, thresh);
		
		System.out.println(" edgels size:       " + edgels.size());
		System.out.println(" edgels float size: " + edgelsf.size());
		
		for(int i=0; i<20; i++){
			System.out.println("d " + edgels.get(i) );
//			System.out.println("f " + edgelsf.get(i) );
			System.out.println(" ");
		}
		
	}
	
	public static void main(String[] args) {

		tryMatching();
		
//		tryMatchRegistration();
		
//		stupdidRadiusSearchTest();
//		
//		edgelTypeValidate();
		
		System.out.println("crack correction finished");
		System.exit(0);
		
	}

}
