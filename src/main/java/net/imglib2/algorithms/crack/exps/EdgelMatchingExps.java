package net.imglib2.algorithms.crack.exps;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import edu.jhu.ece.iacl.utility.ArrayUtil;
import ij.IJ;
import net.imglib2.RandomAccess;
import net.imglib2.RealPoint;
import net.imglib2.algorithm.edge.Edgel;
//import net.imglib2.algorithm.edge.EdgelFloat;
import net.imglib2.algorithm.edge.SubpixelEdgelDetection;
//import net.imglib2.algorithm.edge.SubpixelEdgelDetectionFloat;
import net.imglib2.algorithms.crack.CrackCorrection;
import net.imglib2.algorithms.crack.EdgelMatching;
import net.imglib2.algorithms.crack.EdgelMatching.EdgelPair;
import net.imglib2.collection.KDTree;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.neighborsearch.RadiusNeighborSearchOnKDTree;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;

public class EdgelMatchingExps {

	Img<FloatType> img;
	Img<FloatType> mask;

	CrackCorrection<FloatType> 	cc;
	EdgelMatching<FloatType> 	em;
	
	
	public static <T extends RealType<T> & NativeType<T>>void edgelPatchWrite(CrackCorrection<T> cc, int dsFactor){
		
		String dirBase = "/data-ssd1/john/projects/crackPatching/closeup/patches/matchExpsDepth_"+
				String.format("ds_%d_%d-%d-%d", dsFactor, cc.getPatchSize()[0],cc.getPatchSize()[1],cc.getPatchSize()[2] );
		
		File dirbasef = new File(dirBase);
		if( ! dirbasef.exists() ){ dirbasef.mkdir(); } 
		
		System.out.println(dirBase);
		
//		double[] testPt = new double[]{ 159,253,73 };
		double[] testPt = new double[]{ 173,205,408 }; 
		double[] resPt  = new double[]{ 170, 312, 136 };
		
		ArrayUtil.divide(testPt, dsFactor);
		ArrayUtil.divide(resPt, dsFactor);
		
		double[] diff   = ArrayUtil.clone(resPt);
		ArrayUtil.subtractInPlace(diff, testPt);
		double dist = Math.sqrt(ArrayUtil.sumSquares(diff));
		System.out.println( " dist " + dist);
		
		
		int ei = cc.edgelIdxNearest(testPt);
		Edgel e = cc.getEdgels().get(ei);
		
		String dir = dirBase + File.separator + "test_" + cc.getEdgels().indexOf(e);
		
		File dirf = new File(dir);
		if( ! dirf.exists() ){ dirf.mkdir(); } 
		
//		cc.writeEdgelPatchAndMatches(e, dir);
		cc.writeEdgelDepthPatchAndMatches(e, dir);
		
	}
	
	public static <T extends RealType<T> & NativeType<T>> void testEdgelMatchesWrite(CrackCorrection<T> cc, int downSampleFactor)
	{

//		double[] testPt = new double[]{ 159,253,73 }; 
		double[] testPt = new double[]{ 173,205,408 }; 
		
		ArrayUtil.divide(testPt, downSampleFactor);
		
		int eidx = cc.edgelIdxNearest(testPt);
		Edgel testEdgel = cc.getEdgels().get(eidx);
		
		String edgelCsvFn = "";
		switch ( cc.edgelMatcher.getSearchType() )
		{
		case COUNT:
			edgelCsvFn = "/data-ssd1/john/projects/crackPatching/closeup/matches_ds4/"
					+ "edgelMatches" + cc.edgelMatcher.getEdgelSearchCount() +"_Near_"+ eidx + "_filt.csv";
			break;
		case RADIUS:
			edgelCsvFn = "/data-ssd1/john/projects/crackPatching/closeup/matches_ds4/"
					+ "edgelMatches_Rad-"+ cc.edgelMatcher.getEdgelSearchRadius() + "_Near_"+ eidx + "_filt.csv";
			break;
		default:
			edgelCsvFn = "/data-ssd1/john/projects/crackPatching/closeup/matches_ds4/"
					+ "edgelMatchesNear_"+ eidx + "_filt.csv";
			break;
		}
		
		System.out.println("testEdgel: " + testEdgel );
		try {
			cc.edgelMatcher.computeAndWriteEdgelMatches(edgelCsvFn, testEdgel);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public static void edgelFeatures()
	{
		int downSampleFactor = 4;
		double searchRadius = 100;
//		double searchRadius = 50;
		
		String imgfn = "/data-ssd1/john/projects/crackSegmentation/groundTruth/closeup/img_ds"+downSampleFactor+".tif";
		String maskfn = "/data-ssd1/john/projects/crackSegmentation/groundTruth/closeup/labels_interp_smooth_ds"+downSampleFactor+".tif";
		
//		String featurefn = "/data-ssd1/john/projects/crackPatching/closeup/ds_4_patch_31-31-19/edgelFeatures.csv";
//		int[] patchSize = new int[] { 31, 31, 19 };
		
		String featurefn = "/data-ssd1/john/projects/crackPatching/closeup/ds_4_patch_15-15-9/edgelFeatures.csv";
		int[] patchSize = new int[] { 15, 15, 9 };
		
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
		
		
		EdgelMatching<FloatType> em = new EdgelMatching<FloatType>(img, mask, patchSize);
		em.debugDir = "/data-ssd1/john/projects/crackPatching/cropEdgelMatch";
		em.debugSuffix = "ds"+downSampleFactor;
//		em.debugSuffix = "";
		
//		em.setSearchType( EdgelMatching.SearchTypes.COUNT );
//		em.setEdgelSearchCount(3000);
		
		em.setEdgelSearchRadius( searchRadius / (downSampleFactor) );
		em.setEdgels(cc.getEdgels());
		
		try {
			em.computeAndWriteEdgelsAndFeatures(featurefn);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public static void tryMatchRegistration() 
	{
		String imgfn = "/data-ssd1/john/projects/crackSegmentation/crackVolDown_cp.tif";
		String maskfn = "/data-ssd1/john/projects/crackSegmentation/Labels_ds_interp_cp_smooth.tif";
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
	
		cc.debugOutDir = "/data-ssd1/john/projects/crackPatching/edgelRegPatch";
//		cc.registerEdgelsOrient( cc.getEdgels().get(i), cc.getEdgels().get(j), i);
		
	}
	
	public static void tryMatching() 
	{
//		String imgfn = "/data-ssd1/john/projects/crackSegmentation/crackVolDown_cp.tif";
//		String maskfn = "/data-ssd1/john/projects/crackSegmentation/Labels_ds_interp_cp_smooth.tif";
		
		int downSampleFactor = 4;
		double searchRadius = 100;
//		double searchRadius = 50;
		
//		String imgfn = "/data-ssd1/john/projects/crackSegmentation/groundTruth/closeup/img.tif";
//		String maskfn = "/data-ssd1/john/projects/crackSegmentation/groundTruth/closeup/labels_interp_smooth.tif";
		
		String imgfn = "/data-ssd1/john/projects/crackSegmentation/groundTruth/closeup/img_ds"+downSampleFactor+".tif";
		String maskfn = "/data-ssd1/john/projects/crackSegmentation/groundTruth/closeup/labels_interp_smooth_ds"+downSampleFactor+".tif";
		int[] patchSize = new int[] { 31, 31, 21 };
		
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
		em.debugDir = "/data-ssd1/john/projects/crackPatching/cropEdgelMatch";
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
		
		String maskfn = "/data-ssd1/john/projects/crackSegmentation/Labels_ds_interp_cp_smooth.tif";

		
		ArrayImgFactory<FloatType> ffactory = new ArrayImgFactory<FloatType>();
		ArrayImgFactory<DoubleType> dfactory = new ArrayImgFactory<DoubleType>();
		
		Img<FloatType> mask = ImagePlusAdapter.convertFloat( IJ.openImage(maskfn) );
	
		float thresh = 40.0f;
		
		ArrayList<Edgel> 		edgels  = SubpixelEdgelDetection.getEdgels(mask, ffactory, thresh);
//		ArrayList<EdgelFloat> 	edgelsf = SubpixelEdgelDetectionFloat.getEdgels(mask, ffactory, thresh);
		
		System.out.println(" edgels size:       " + edgels.size());
//		System.out.println(" edgels float size: " + edgelsf.size());
		
		for(int i=0; i<20; i++){
			System.out.println("d " + edgels.get(i) );
//			System.out.println("f " + edgelsf.get(i) );
			System.out.println(" ");
		}
		
	}
	
	
	public static <T extends RealType<T> & NativeType<T>> void makeLapEdgeImg( CrackCorrection<T> cc, int dsFactor){
		
		ArrayList<Edgel> eList = cc.getEdgels();
		
		Img<T> edgeImg = cc.getImg().factory().create(cc.getImg(), cc.getImg().firstElement());
		RandomAccess<T> ra = edgeImg.randomAccess();
		
		double[] pos = new double[ cc.getImg().numDimensions() ];
		int[] posi = null;
		
		int i = 0;
		for( Edgel e : eList )
		{
			double depth = cc.computeDepthLap( e );
			e.localize( pos );

			double[] grd = ArrayUtil.clone(e.getGradient());
			ArrayUtil.normalizeLengthInPlace(grd);
			ArrayUtil.multiply(grd, depth);
			
			ArrayUtil.addInPlace( pos, grd );
		
			posi = ArrayUtil.toIntRound( pos );
			ra.setPosition( posi );
			ra.get().setOne();
			
			i++;
			if( i % 10 == 0 ){
				System.out.println(" edgel " + i + " of " + eList.size() );
			}
		}
		
		ImgOps.writeFloat( edgeImg, "/data-ssd1/john/projects/crackPatching/closeup/lapEdgeImg.tif" );
	}
	
	
	public static void main(String[] args) {
		
		int downSampleFactor = 4;
		double searchRadius = 100;
		int    searchCount  = 2000;
		
		String imgfn = "/data-ssd1/john/projects/crackSegmentation/groundTruth/closeup/img_ds"+downSampleFactor+".tif";
		String maskfn = "/data-ssd1/john/projects/crackSegmentation/groundTruth/closeup/labels_interp_smooth_ds"+downSampleFactor+".tif";
		Img<FloatType> img =  ImagePlusAdapter.convertFloat( IJ.openImage(imgfn) );
		Img<FloatType> mask = ImagePlusAdapter.convertFloat( IJ.openImage(maskfn) );
		
//		int[] patchSize = new int[] { 45, 45, 19 };
		int[] patchSize = new int[] { 31, 31, 19 };
//		int[] patchSize = new int[] { 15, 15, 13 };
		
		CrackCorrection<FloatType> cc = new CrackCorrection<FloatType>(
				img, mask, patchSize);
		
		EdgelMatching<FloatType> em = new EdgelMatching<FloatType>(img, mask, patchSize);
		em.debugDir = "/data-ssd1/john/projects/crackPatching/cropEdgelMatch";
		em.debugSuffix = "ds"+downSampleFactor;
		
		cc.edgelMatcher = em;
		
		cc.edgelMatcher.setEdgelSearchRadius( searchRadius / (downSampleFactor) );
		
//		cc.edgelMatcher.setSearchType( EdgelMatching.SearchTypes.COUNT );
//		cc.edgelMatcher.setEdgelSearchCount( searchCount );
		
		cc.computeEdgels();
		cc.edgelMatcher.setEdgels(cc.getEdgels());
		
//		makeLapEdgeImg( cc, downSampleFactor );
		
		testEdgelMatchesWrite( cc, downSampleFactor );
		
//		edgelPatchWrite( cc, downSampleFactor );
		
//		edgelFeatures();
		
//		tryMatching();
		
//		tryMatchRegistration();
		
//		stupdidRadiusSearchTest();
//		
//		edgelTypeValidate();
		
		System.out.println("execution finished");
		System.exit(0);
		
	}

}
