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

public class EdgelMatchingExps {

	public static void tryMatching() 
	{
		String imgfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/crackVolDown_cp.tif";
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_smooth.tif";
		int[] patchSize = new int[] { 19, 19, 13 };
		
		Img<FloatType> img =  ImagePlusAdapter.convertFloat( IJ.openImage(imgfn) );
		Img<FloatType> mask = ImagePlusAdapter.convertFloat( IJ.openImage(maskfn) );

		CrackCorrection<FloatType, FloatType> cc = new CrackCorrection<FloatType, FloatType>(
				img, mask, patchSize);
		
		cc.computeEdgels();
		
//		ArrayList<Edgel> edgels = cc.getEdgels();
//		int i =0;
//		for ( Edgel e : edgels ){
//			System.out.println( ArrayUtil.printArray( e.getPosition() ));
//			
//			if(i>50) break;
//		}
		
		EdgelMatching<FloatType> em = new EdgelMatching<FloatType>(img, mask, patchSize);
		em.setEdgels(cc.getEdgels());
		
//		em.setEdgelSearchRadius(20);
		
		em.setSearchType( EdgelMatching.SearchTypes.COUNT );
		em.setEdgelSearchCount(50);
		
		em.computeAllAffinities();
		
		
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
			System.out.println("d " + ArrayUtil.printArray( edgels.get(i).getPosition()));
			System.out.println("f " + ArrayUtil.printArray( edgelsf.get(i).getPosition()));
			System.out.println(" ");
		}
		
	}
	
	public static void main(String[] args) {

		tryMatching();
		
//		stupdidRadiusSearchTest();
//		
//		edgelTypeValidate();
		
		System.out.println("crack correction finished");
		System.exit(0);
		
	}

}
