package net.imglib2.algorithms.crack;

import ij.IJ;
import ij.ImagePlus;

import java.util.ArrayList;
import java.util.Collections;

import net.imglib2.io.*;
import net.imglib2.type.*;
import net.imglib2.type.numeric.integer.*;
import net.imglib2.type.numeric.real.*;
import edu.jhu.ece.iacl.algorithms.MGDM.MgdmDecomposition;
import edu.jhu.ece.iacl.utility.ArrayUtil;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.util.ImgUtil;

public class CrackToyExamples {

	public static void crackToy1()
	{
		String crackMaskFn = "/groups/jain/home/bogovicj/projects/crackSegmentation/toy/toyCrack.tif";
		String depthFn = "/groups/jain/home/bogovicj/projects/crackSegmentation/toy/toyCrack_depth.tif";
//		String patchFn = "/groups/jain/home/bogovicj/projects/crackSegmentation/toy/toyCrack_patch.tif";
		String patch2Fn = "/groups/jain/home/bogovicj/projects/crackSegmentation/toy/toyCrack_patch2.tif";

		IntType type = new IntType();
		ArrayImgFactory< IntType > ifactory = new ArrayImgFactory< IntType >();
		ArrayImgFactory< FloatType > ffactory = new ArrayImgFactory< FloatType >();
		
		Img<IntType> mask = null;
		Img<FloatType> img  = null;
		try 
		{
			mask = new ImgOpener().openImg( crackMaskFn , ifactory, type );
			img = ImgUtil.createGradientImgY((int)mask.dimension(0), 
					(int)mask.dimension(1), (int)mask.dimension(2), new FloatType());
		}
		catch (ImgIOException e)
		{
			e.printStackTrace();
		}
		
		ArrayList<Integer> uniqueList = ImgUtil.uniqueInt(mask);
		Collections.sort(uniqueList);
		System.out.println("uniqueList: " + uniqueList);

		CrackCorrection<FloatType, IntType> cc = 
				new CrackCorrection<FloatType, IntType>( 
						img, mask );

		float[] tgtPos  = new float[]{9f, 11f, 3f};
		int[] patchSize = new int[]{ 7, 7 };
		
		cc.computeEdgels();
		int i = cc.edgelIdxNearest( tgtPos );
		Edgel edgel = cc.edgels.get(i);
	
	
//		System.out.println(" \nedgel pos : " + ArrayUtil.printArray(edgel.getPosition()));
//		System.out.println(" edgel grad: " + ArrayUtil.printArray(edgel.getGradient()));
//		System.out.println(" edgel mag : " + edgel.getMagnitude());
//		
//		Img<FloatType> depthPatch = cc.computeLocalCrackDepthDistFunc(edgel, patchSize);
//		Img<FloatType> imgPatch   = cc.edgelToImage(edgel, img, patchSize);
//		
//		write( depthPatch, depthFn );
//		write( imgPatch, patchFn );
		
		Img<FloatType> imgPatch = cc.computeCrackDepthNormal(edgel, patchSize);
//		write( imgPatch, patch2Fn );
		
	}
	

	public static <T extends NativeType<T>> void write(Img<T> img, String fn)
	{
		try
		{
			ImagePlus ipdp = ImgUtil.toImagePlus( img );
			IJ.save(ipdp, fn);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	
	public static void testMgdm(){
		
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/toy/toyCrack.tif";
		
		IntType type = new IntType();
		ArrayImgFactory< IntType > factory = new ArrayImgFactory< IntType >();
		
		Img<IntType> mask = null;
		try 
		{
			mask = new ImgOpener().openImg( maskfn , factory, type );
		}
		catch (ImgIOException e)
		{
			e.printStackTrace();
		}
		
		ImgUtil.printNumNonZero(mask);
		ArrayList<Integer> uniqueList = ImgUtil.uniqueInt(mask);
		
		System.out.println( "uniqueList: " + uniqueList );
		
		int[][][] maskInt = ImgUtil.toIntArray3d(mask);
		
//		System.out.println("maskInt:\n"+ maskInt.length + " " + maskInt[0].length + " " + maskInt[0][0].length);
//		System.out.println("maskInt:\n"+ ArrayUtil.printArray(maskInt));
		
		boolean[][][] m = new boolean[maskInt.length][maskInt[0].length][maskInt[0][0].length];
		ArrayUtil.fill(m, true);
		

		MgdmDecomposition mgdm = new MgdmDecomposition(maskInt, 2, m);
//		float[][][] d1 = mgdm.exportDistanceForObject3d(1);
//		System.out.println("dist:\n"+ ArrayUtil.printArray(d1));
		
		float[][][][] dists = mgdm.exportDistances3d();
		System.out.println("dists:\n" );
		print(dists);
		
//		int x = 9;
//		int y = 11;
//		int z = 3;
//		float dist0 = mgdm.exportDistanceForObject3d(x, y, z, 0);
//		float dist1 = mgdm.exportDistanceForObject3d(x, y, z, 1);
//
//		System.out.println("dist:\n" + dist0);
//		System.out.println("dist:\n" + dist1);
		
		
//		Img<FloatType> distXfmImg = null;
//
//		logger.info( "writing patches" );
//		ArrayImgFactory<FloatType> factory = new ArrayImgFactory<FloatType>();
//
//		distXfmImg = factory.create(img, new FloatType(img
//				.firstElement().getRealFloat()));
//		ImgUtil.copyToImg(distXfmImg, d1);
		
		
	}
	
//	public static void testDistanceMap(){
//		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/toy/toyCrack.tif";
//		
//		IntType type = new IntType();
//		ArrayImgFactory< IntType > factory = new ArrayImgFactory< IntType >();
//		
//		Img<IntType> mask = null;
//		try 
//		{
//			mask = new ImgOpener().openImg( maskfn , factory, type );
//		}
//		catch (ImgIOException e)
//		{
//			e.printStackTrace();
//		}
//		
//		DistanceMap dm;
//	}
	
	public static void print( float[][][][] in ){

		String out = "";
		for(int n=0; n<in[0][0][0].length; n++){
			for(int k=0; k<in[0][0].length; k++){

				for(int i=0; i<in.length; i++){
					for(int j=0; j<in[0].length; j++){
						out += in[i][j][k][n] +" ";
					}
					out +="\n";
				}
				out += "****************\n";
			}
			out += "****************\n";
		}

		System.out.println(out);
	}
	
	public static void main(String[] args) {

//		crackToy1();
		testMgdm();
		
		System.out.println("done");
		System.exit(0);
	}

}
