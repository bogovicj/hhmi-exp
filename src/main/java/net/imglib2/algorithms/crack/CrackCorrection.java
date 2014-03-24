package net.imglib2.algorithms.crack;

import java.util.ArrayList;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.ejml.data.DenseMatrix64F;

import edu.jhu.ece.iacl.algorithms.MGDM.MgdmDecomposition;
import edu.jhu.ece.iacl.utility.ArrayUtil;
import ij.IJ;
import ij.ImagePlus;
import net.imglib2.*;
import net.imglib2.img.imageplus.*;
import net.imglib2.algorithm.edge.*;
import net.imglib2.algorithms.moments.ImgMoment;
import net.imglib2.img.*;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.*;
import net.imglib2.type.numeric.integer.*;
import net.imglib2.type.numeric.real.*;
//import net.imglib2.algorithm.integral.*;

public class CrackCorrection<T extends RealType<T>, B extends AbstractIntegerType<B>> {

	Img<T> img;
	Img<B> mask;

	//IntegralImgDouble<B> maskIntImg;
	ArrayList<Edgel> edgels;
	
	Img<B> imgChunks;
	
	protected static Logger logger = LogManager.getLogger(CrackCorrection.class.getName());

	public CrackCorrection (){ }

	public CrackCorrection( Img<T> img, Img<B> mask ){
		this.img = img;
		this.mask = mask;
	}

	public void genImgChunks(){

		if(imgChunks == null){
			ImgFactory<B> factory = mask.factory();	
			imgChunks = factory.create(mask, mask.firstElement());
		}
		//DenseMatrix64F mtx;
	}

	public void computeEdgels(){
		
		System.out.println("mask " + mask);
		edgels = SubpixelEdgelDetection.getEdgels(mask, mask.factory(), 2f);
		System.out.println("num edgels " + edgels.size());
		
		//System.out.println("\nfirst edgel: ");
		//System.out.println(" pos: " + ArrayUtil.printArray( edgels.get(0).getPosition() ));
		//System.out.println(" grd: " + ArrayUtil.printArray( edgels.get(0).getGradient() ));
		//System.out.println(" mag: " + edgels.get(0).getMagnitude() );
		
	}

	public void crackBoundaries(){
		// compute crack orientation:
		ImgMoment imgmom = new ImgMoment();
		double[] or = imgmom.orientation( mask );	
		imgmom.orientationEvalsEvecs( or );

		double[] cent = imgmom.getM1();
		ArrayUtil.divide(cent, imgmom.getM0());

		//System.out.println(" pos: " + ArrayUtil.printArray( ArrayUtil.reshape2D( or , 3, 3, false )));
		System.out.println(" cent:\n " + ArrayUtil.printArray( cent ));
		System.out.println(" orEvecs:\n " + ArrayUtil.printArray( ArrayUtil.reshape2D( imgmom.getEvecs(), 3, 3, false )));	


	}

	public static <T extends RealType<T>> int[][][] toIntArray3d(Img<T> img){
		int[][][] out = new int[(int)img.dimension(0)][(int)img.dimension(1)][(int)img.dimension(2)];
		Cursor<T> cursor = img.localizingCursor();
		int[] pos = new int[3];
		while(cursor.hasNext()){
			cursor.next();
			cursor.localize(pos);
			out[pos[0]][pos[1]][pos[2]] = (int)(cursor.get().getRealDouble());
		}
		return out;
	}
	
	public static <T extends RealType<T>> boolean[][][] toBooleanArray3dNeg(Img<T> img){
		boolean[][][] out = new boolean[(int)img.dimension(0)][(int)img.dimension(1)][(int)img.dimension(2)];
		Cursor<T> cursor = img.localizingCursor();
		int[] pos = new int[3];
//		long num = 0;
		while(cursor.hasNext()){
			cursor.next();
			cursor.localize(pos);
			out[pos[0]][pos[1]][pos[2]] = (cursor.get().getRealDouble() < 0.5);
			
//			if(out[pos[0]][pos[1]][pos[2]]){
//				num++;
//			}
		}
		
//		System.out.println("num: " + num);
		
		return out;
	}

	public static <T extends RealType<T>> float[][][] toFloatArray3d(Img<T> img){
		float[][][] out = new float[(int)img.dimension(0)][(int)img.dimension(1)][(int)img.dimension(2)];
		Cursor<T> cursor = img.localizingCursor();
		int[] pos = new int[3];
		while(cursor.hasNext()){
			cursor.next();
			cursor.localize(pos);
			out[pos[0]][pos[1]][pos[2]] = (cursor.get().getRealFloat());
		}
		return out;
	}
	
	public static <T extends RealType<T>> void copyToImg4d(Img<T> img, int[][][][] in){
		
		Cursor<T> cursor = img.localizingCursor();
		int[] pos = new int[4];
		
		while(cursor.hasNext()){

			cursor.next();
			cursor.localize(pos);
			cursor.get().setReal( (double) in[pos[0]][pos[1]][pos[2]][pos[3]] );

		}
	}
	
	public static <T extends RealType<T>> void copyToImg4d(Img<T> img, float[][][][] in){
		
		Cursor<T> cursor = img.localizingCursor();
		int[] pos = new int[4];
		
		while(cursor.hasNext()){

			cursor.next();
			cursor.localize(pos);
			cursor.get().setReal(  in[pos[0]][pos[1]][pos[2]][pos[3]] );

		}
	}

	public void computeLocalCrackDepth(){
		logger.debug("computeLocalCrackDepth");	

	}
	
	public static < T extends NativeType< T >> ImagePlusImg<T, ?> copyToImagePlus(Img<T> img) throws Exception {


		ImagePlusImgFactory<T> factory = new ImagePlusImgFactory<T>();
		ImagePlusImg<T, ?> ipImg = factory.create(img, img.firstElement());

		System.out.println("create image plus of type: " + img.firstElement().getClass());
		System.out.println("result is of type: " + ipImg.firstElement().getClass());


		Cursor<T> c_in  = img.cursor();
		RandomAccess<T> ra = ipImg.randomAccess();

		while(c_in.hasNext()){
			c_in.fwd();
			ra.setPosition(c_in);
			ra.get().set(c_in.get());

		}

		return ipImg;

	}

	public static < T extends NativeType< T >> ImagePlus toImagePlus(Img<T> img) throws Exception{                 
		return copyToImagePlus(img).getImagePlus();
	} 

	   
	public static void main(String[] args){

		String imgfn  = "/groups/jain/home/bogovicj/projects/crackSegmentation/crackVolDown_cp.tif";
		String maskfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp.tif";
		String mgdmfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/LabelsMgdm_ds_interp_cp_v2.tif";
		
		String maskWSfn = "/groups/jain/home/bogovicj/projects/crackSegmentation/Labels_ds_interp_cp_manPaint2.tif";

		
		ImagePlus imgip = IJ.openImage(imgfn);
		ImagePlus maskip = IJ.openImage(maskfn);

//		System.out.println("imgip " + imgip);
//		System.out.println("maskip" + maskip);

		Img<FloatType> img = ImagePlusAdapter.convertFloat(imgip);
		ByteImagePlus<UnsignedByteType> mask = ImagePlusAdapter.wrapByte(maskip); 

//		System.out.println("img " + img);
//		System.out.println("mask " + mask);

		CrackCorrection<FloatType, UnsignedByteType> cc = 
				new CrackCorrection<FloatType, UnsignedByteType>( 
						img, mask );

		//cc.computeEdgels();
		//cc.crackBoundaries();

		boolean[][][]   mskArray = toBooleanArray3dNeg(mask);
		//float[][][] imgArray = toFloatArray3d(img);
		

	      System.out.println("nnz in mask: \n" +
	               ArrayUtil.nnz(ArrayUtil.reshape1D(mskArray, true))
	            );
	      
	      int[][][]  sideArray = new int[mskArray.length][mskArray[0].length][mskArray[0][0].length];
	      
	      sideArray[0][0][13]     = 1;
	      sideArray[287][437][13] = 2;

	      MgdmDecomposition mgdm = new MgdmDecomposition(sideArray,3, mskArray);
	      
	      int[][][][] mgdmLabels = mgdm.exportAllLabels3d();
	      
	      Img<UnsignedByteType> mgdmLImg = mask.factory().create(
	            new long[]{ mask.dimension(0), mask.dimension(1),
	                     mask.dimension(2), 3 },
	            new UnsignedByteType());
	      
//	    Img<FloatType> d1img = img.factory().create(img, img.firstElement());
//	    Img<FloatType> d2img = img.factory().create(img, img.firstElement());
//	   
//	    float[][][] d1 = mgdm.exportDistanceForObject3d(1);
//	    float[][][] d2 = mgdm.exportDistanceForObject3d(2);
	      
	      try {
	         
	         copyToImg4d(mgdmLImg, mgdmLabels);
	         ImagePlus ipMgdm = toImagePlus(mgdmLImg);
	         IJ.save(ipMgdm, mgdmfn);
	         
//	       copyToImg4d(d1img, d1);
//	       ImagePlus ipMgdm = toImagePlus(d1img);
//	       IJ.save(ipMgdm, d1fn);
//	       
//	       copyToImg4d(d2img, d2);
//	       ipMgdm = toImagePlus(d2img);
//	       IJ.save(ipMgdm, d2fn);
	         
	      } catch (Exception e) {
	         e.printStackTrace();
	      }


		System.out.println("crack correction finished");
		System.exit(0);
	}


}
