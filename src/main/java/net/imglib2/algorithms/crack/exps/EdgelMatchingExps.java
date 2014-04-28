package net.imglib2.algorithms.crack.exps;

import ij.IJ;
import net.imglib2.algorithms.crack.CrackCorrection;
import net.imglib2.algorithms.crack.EdgelMatching;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
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
		
		EdgelMatching<FloatType> em = new EdgelMatching<FloatType>(img, mask, patchSize);
		em.setEdgels(cc.getEdgels());
		
		em.computeAllAffinities();
		
		
	}
	
	
	public static void main(String[] args) {

		tryMatching();
		
		
		System.out.println("crack correction finished");
		System.exit(0);
		
	}

}
