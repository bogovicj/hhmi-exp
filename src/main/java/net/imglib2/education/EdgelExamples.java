package net.imglib2.education;

import java.util.ArrayList;

import edu.jhu.ece.iacl.utility.ArrayUtil;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithm.edge.SubpixelEdgelDetection;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgUtil;

public class EdgelExamples {

	public static Img<FloatType> genImage()
	{
		 FloatType t = new FloatType();
		 Img<FloatType> img = ImgUtil.createEdgeImg( new int[]{64,64,12}, new double[]{1,2,0}, t, 2);
		 
		 String fnOut = "/Users/bogovicj/Documents/projects/crackStitching/edgeTestImg.tif";
		 ImgUtil.write(img, fnOut);
		 
		 return img;
	}
	
	public static void main(String[] args) {
		System.out.println("Starting");
		
		Img<FloatType> img = genImage();
		float minMag = 0.1f;
		
		ArrayList<Edgel> edgels = SubpixelEdgelDetection.getEdgels( img, img.factory(), minMag);
		System.out.println("found " + edgels.size() + " edgels");
		
//		Edgel edgel = edgels.get(0);
//		System.out.println("edgel: " + edgel);
//		
//		System.out.println(" edgel pos : " + ArrayUtil.printArray(edgel.getPosition()));
//		System.out.println(" edgel grad: " + ArrayUtil.printArray(edgel.getGradient()));
//		System.out.println(" edgel mag : " + edgel.getMagnitude());
		
		System.out.println("Done");
		System.exit(0);
	}

}
