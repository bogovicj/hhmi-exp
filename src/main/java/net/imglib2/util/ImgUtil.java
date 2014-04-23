package net.imglib2.util;

import java.util.ArrayList;
import java.util.HashSet;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import edu.jhu.ece.iacl.utility.ArrayUtil;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.imageplus.ImagePlusImg;
import net.imglib2.img.imageplus.ImagePlusImgFactory;
import net.imglib2.ops.operation.randomaccessibleinterval.unary.DistanceMap;
import net.imglib2.type.NativeType;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.AbstractIntegerType;
import net.imglib2.type.numeric.integer.ByteType;
import net.imglib2.type.numeric.integer.GenericByteType;
import net.imglib2.type.numeric.integer.GenericShortType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.integer.ShortType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;

public class ImgUtil {

	public static Logger logger = LogManager.getLogger(ImgUtil.class.getName());
	
	public static <T extends NativeType<T>> void write(Img<T> img, String fn)
	{
		try
		{
			ImagePlus ipdp = ImgUtil.copyToImagePlus(img).getImagePlus();
			IJ.save(ipdp, fn);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	
	public static void writeFloat(Img<FloatType> img, String fn)
	{
		try
		{
			ImagePlus ipdp = mimicImagePlusFloat( img, "name" );
			
			if(img.numDimensions()==2){
				copyToImageProcessor2dFloat(img, ipdp.getProcessor());
			}else if(img.numDimensions()==3){
				copyToImageProcessor3dFloat(img, ipdp);
			}
			
			IJ.save(ipdp, fn);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	public static <T extends GenericByteType<T>> void writeByte(Img<T> img, String fn)
	{
		try
		{
			ImagePlus ipdp = mimicImagePlusByte( img, "name" );
			if(img.numDimensions()==2){
				copyToImageProcessor2dByte(img, ipdp.getProcessor());
			}else if(img.numDimensions()==3){
				copyToImageProcessor3dByte(img, ipdp);
			}
			IJ.save(ipdp, fn);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	public static <T extends GenericShortType<T>> void writeShort(Img<T> img, String fn)
	{
		try
		{
			ImagePlus ipdp = mimicImagePlusShort( img, "name" );
			if(img.numDimensions()==2){
				copyToImageProcessor2dShort(img, ipdp.getProcessor());
			}else if(img.numDimensions()==3){
				copyToImageProcessor3dShort(img, ipdp);
			}
			IJ.save(ipdp, fn);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	public static ImagePlus mimicImagePlusFloat(Img<FloatType> img, String name)
	{
		ImagePlus ip = null;
		if(img.numDimensions()==2){
			FloatProcessor fp = new FloatProcessor(
					(int)img.dimension(0),(int)img.dimension(1));
			ip = new ImagePlus(name, fp);
			
		}else if(img.numDimensions()==3){
			
			ImageStack stack = ImageStack.create(
					(int)img.dimension(0),
					(int)img.dimension(1),
					(int)img.dimension(2),
					32);
			ip = new ImagePlus(name, stack);
		}
		
		return ip;
	}
	public static <T extends GenericByteType<T>> ImagePlus mimicImagePlusByte(
			Img<T> img, String name)
	{
		ImagePlus ip = null;
		if(img.numDimensions()==2){
			ByteProcessor fp = new ByteProcessor(
					(int)img.dimension(0),(int)img.dimension(1));
			ip = new ImagePlus(name, fp);
			
		}else if(img.numDimensions()==3){
			
			ImageStack stack = ImageStack.create(
					(int)img.dimension(0),
					(int)img.dimension(1),
					(int)img.dimension(2),
					8);
			ip = new ImagePlus(name, stack);
		}
		
		return ip;
	}
	public static <T extends GenericShortType<T>> ImagePlus mimicImagePlusShort(
			Img<T> img, String name)
	{
		ImagePlus ip = null;
		if(img.numDimensions()==2){
			ByteProcessor fp = new ByteProcessor(
					(int)img.dimension(0),(int)img.dimension(1));
			ip = new ImagePlus(name, fp);
			
		}else if(img.numDimensions()==3){
			
			ImageStack stack = ImageStack.create(
					(int)img.dimension(0),
					(int)img.dimension(1),
					(int)img.dimension(2),
					16);
			ip = new ImagePlus(name, stack);
		}
		
		return ip;
	}
	
	public static <T extends GenericByteType<T>> void copyToImageProcessor2dByte(Img<T> img, ImageProcessor ip) throws Exception {
		
		Cursor<T> c_in  = img.cursor();
		
		while(c_in.hasNext()){
			c_in.fwd();
			
			ip.set(c_in.getIntPosition(0), c_in.getIntPosition(1), 
					c_in.get().getInteger());
			
		}
	}
	public static <T extends GenericShortType<T>> void copyToImageProcessor2dShort(Img<T> img, ImageProcessor ip) throws Exception {
		
		Cursor<T> c_in  = img.cursor();
		
		while(c_in.hasNext()){
			c_in.fwd();
			
			ip.set(c_in.getIntPosition(0), c_in.getIntPosition(1), 
					c_in.get().getInteger());
			
		}
	}
	public static void copyToImageProcessor2dFloat(Img<FloatType> img, ImageProcessor ip) throws Exception {
		
		Cursor<FloatType> c_in  = img.cursor();
		
		while(c_in.hasNext()){
			c_in.fwd();
			
			ip.setf(c_in.getIntPosition(0), c_in.getIntPosition(1), 
					c_in.get().getRealFloat());
			
		}
	}
	
	public static void copyToImageProcessor3dFloat(Img<FloatType> img, ImagePlus ip) throws Exception {
		
		Cursor<FloatType> c_in  = img.cursor();
		
		while(c_in.hasNext()){
			c_in.fwd();
			
			ImageProcessor fp = ip.getStack().getProcessor(c_in.getIntPosition(2)+1);
			
			fp.setf(c_in.getIntPosition(0), c_in.getIntPosition(1), 
					c_in.get().getRealFloat());
			
		}
	}
	public static <T extends GenericByteType<T>> void copyToImageProcessor3dByte(Img<T> img, ImagePlus ip) throws Exception {
		
		Cursor<T> c_in  = img.cursor();
		
		while(c_in.hasNext()){
			c_in.fwd();
			
			ImageProcessor fp = ip.getStack().getProcessor(c_in.getIntPosition(2)+1);
			
			fp.set(c_in.getIntPosition(0), c_in.getIntPosition(1), 
					c_in.get().getInteger());
			
		}
	}
	public static <T extends GenericShortType<T>> void copyToImageProcessor3dShort(Img<T> img, ImagePlus ip) throws Exception {
		
		Cursor<T> c_in  = img.cursor();
		
		while(c_in.hasNext()){
			c_in.fwd();
			
			ImageProcessor fp = ip.getStack().getProcessor(c_in.getIntPosition(2)+1);
			
			fp.set(c_in.getIntPosition(0), c_in.getIntPosition(1), 
					c_in.get().getInteger());
			
		}
	}
	public static <B extends AbstractIntegerType<B>> Img<FloatType> signedDistance(Img<B> mask){
		
		DistanceMap<B> dm = new DistanceMap<B>();
		ArrayImgFactory<FloatType> ffactory = new ArrayImgFactory<FloatType>();
		Img<FloatType> sdf = ffactory.create(mask, new FloatType());
		Img<FloatType> sdfi = ffactory.create(mask, new FloatType());
		
		
//		ImgUtil.numNonZero(mask);
		
		// inside region
		dm.compute(mask, sdf);
		
		// negate mask
		Img<B> maskInv = ImgUtil.thresholdMap(mask, 1, false);
//		ImgUtil.numNonZero(maskInv);
		
		// oustide region
		dm.compute(maskInv, sdfi);
		
		 Cursor<FloatType> cursor = sdf.cursor();
		 Cursor<FloatType> cursori = sdfi.cursor();
		 while(cursor.hasNext()){
			 cursor. fwd();
			 cursori.fwd();
			 
			 // Distance map returns squared distance
			 // 
			 double inVal  = Math.sqrt(cursor.get().getRealDouble());
			 double outVal = Math.sqrt(cursori.get().getRealDouble());
			 
			 if(inVal > 0){  inVal -= 0.5; }
			 if(outVal > 0){ outVal -= 0.5; }
			 
			 inVal *= -1; 						 // negate inside
			 cursor.get().set( (float)(inVal + outVal) ); // add outside
		 }
		
		return sdf;
	}
	
	public static <T extends Type<T>> void fill(Img<T> img, T value){
		Cursor<T> cursor = img.cursor();
		while(cursor.hasNext()){
			cursor.next().set(value);
		}
	}
	
	public static <S extends RealType<S>> int numNonZero(Img<S> img){
		Cursor<S> c = img.cursor();
		int num = 0;
		while(c.hasNext()){
			S val = c.next();
			if( val.getRealDouble() != 0 ){
				num++;
			}
		}
//		System.out.println(img + " nnz: " + num );
		return num;
	}
	

	public static <S extends RealType<S>> void printCoordNonZero(Img<S> img){
		Cursor<S> c = img.cursor();
		int num = 0;
		int[] pos = new int[img.numDimensions()];
				
		while(c.hasNext()){
			S val = c.next();
			if( val.getRealDouble() != 0 ){
				c.localize(pos);
				System.out.println(" val of : " + val + " at: " + ArrayUtil.printArray(pos));
				num++;
			}
		}
		System.out.println(img + " nnz: " + num );
	}

   public static <T extends RealType<T>> int[][][] toIntArray3d(Img<T> img){
      int[][][] out = new int[(int)img.dimension(0)][(int)img.dimension(1)][(int)img.dimension(2)];
      Cursor<T> cursor = img.localizingCursor();
      int[] pos = new int[3];
      while(cursor.hasNext()){
         cursor.fwd();
         cursor.localize(pos);
         out[pos[0]][pos[1]][pos[2]] = (int)(cursor.get().getRealDouble());
      }
      return out;
   }
   
   public static <T extends RealType<T>> boolean[][][] toBooleanArray3dNeg(Img<T> img){
      boolean[][][] out = new boolean[(int)img.dimension(0)][(int)img.dimension(1)][(int)img.dimension(2)];
      Cursor<T> cursor = img.localizingCursor();
      int[] pos = new int[3];
      while(cursor.hasNext()){
         cursor.next();
         cursor.localize(pos);
         out[pos[0]][pos[1]][pos[2]] = (cursor.get().getRealDouble() < 0.5);
      }
      
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
   
   public static <T extends RealType<T>> void copyToImg(Img<T> img, int[][][] in){
      
      Cursor<T> cursor = img.localizingCursor();
      int[] pos = new int[3];
      while(cursor.hasNext()){
         cursor.next();
         cursor.localize(pos);
         cursor.get().setReal( (double) in[pos[0]][pos[1]][pos[2]] );
      }
   }
   
   public static <T extends RealType<T>> void copyToImg(Img<T> img, float[][][] in){
      
      Cursor<T> cursor = img.localizingCursor();
      int[] pos = new int[3];
      
      while(cursor.hasNext()){
         cursor.next();
         cursor.localize(pos);
         cursor.get().setReal(  in[pos[0]][pos[1]][pos[2]]);
      }
   }
   
   public static <T extends RealType<T>> void copyToImg(Img<T> img, int[][][][] in){
      
      Cursor<T> cursor = img.localizingCursor();
      int[] pos = new int[4];
      while(cursor.hasNext()){
         cursor.next();
         cursor.localize(pos);
         cursor.get().setReal( (double) in[pos[0]][pos[1]][pos[2]][pos[3]] );
      }
   }
   
   public static <T extends RealType<T>> void copyToImg(Img<T> img, float[][][][] in){
      
      Cursor<T> cursor = img.localizingCursor();
      int[] pos = new int[4];
      
      while(cursor.hasNext()){
         cursor.next();
         cursor.localize(pos);
         cursor.get().setReal(  in[pos[0]][pos[1]][pos[2]][pos[3]] );
      }
   }
   

   public static <T extends NativeType<T> & RealType<T>> Img<T> createEdgeImg(int[] sz, double[] w,  T t, double sigma){
      ArrayImgFactory<T> factory = new ArrayImgFactory<T>();
      Img<T> out = factory.create( sz, t);
      
      double[] ctr = ArrayUtil.toDouble(sz);
      ArrayUtil.divide(ctr, 2);
      
      System.out.println(" ctr = " + ArrayUtil.printArray(ctr));
      
      Cursor<T> c = out.localizingCursor();
      double[] pos = new double[3];
      while(c.hasNext()){
         T val = c.next();
         c.localize(pos);
         double[] pt = ArrayUtil.subtract(pos, ctr);
         double[] res = ArrayUtil.multiply( w , pt);
         val.setReal( 
               sigmoid( ArrayUtil.sum(res), sigma )
            );
      }
      
      return out;
   }
   
   public static double sigmoid(double x, double sigma)
   {
	   return (1 / (1 + Math.exp( - sigma * x )));
   }
   
   public static <T extends NativeType<T> & RealType<T>> Img<T> createGradientImg(int[] sz, double[] w, T t){
      ArrayImgFactory<T> factory = new ArrayImgFactory<T>();
      Img<T> out = factory.create( sz, t);
      
      Cursor<T> c = out.localizingCursor();
      int[] pos = new int[3];
      while(c.hasNext()){
         T val = c.next();
         c.localize(pos);
         double[] res = ArrayUtil.multiply( w , ArrayUtil.toDouble(pos));
         val.setReal( 
               ArrayUtil.sum(res)
            );
      }
      
      return out;
   }
   
   public static <T extends NativeType<T> & RealType<T>> Img<T> createGradientImgX(int width, int height, int depth, T t){
      ArrayImgFactory<T> factory = new ArrayImgFactory<T>();
      Img<T> out = factory.create(new int[]{width,height,depth}, t);
      
      Cursor<T> c = out.localizingCursor();
      int[] pos = new int[3];
      while(c.hasNext()){
         T val = c.next();
         c.localize(pos);
         val.setReal(pos[0]);
      }
      
      return out;
   }

   public static <T extends NativeType<T> & RealType<T>> Img<T> createGradientImgY(int width, int height, int depth, T t){
      ArrayImgFactory<T> factory = new ArrayImgFactory<T>();
      Img<T> out = factory.create(new int[]{width,height,depth}, t);
      
      Cursor<T> c = out.localizingCursor();
      int[] pos = new int[3];
      while(c.hasNext()){
         T val = c.next();
         c.localize(pos);
         val.setReal(pos[1]);
      }
      
      return out;
   }
   
   public static <T extends NativeType<T> & RealType<T>> Img<T> createGradientImgZ(int width, int height, int depth, T t){
      ArrayImgFactory<T> factory = new ArrayImgFactory<T>();
      Img<T> out = factory.create(new int[]{width,height,depth}, t);
      
      Cursor<T> c = out.localizingCursor();
      int[] pos = new int[3];
      while(c.hasNext()){
         T val = c.next();
         c.localize(pos);
         val.setReal(pos[2]);
      }
      
      return out;
   }
   
   public static < T extends RealType< T >> Img<T> threshold(Img<T> img, double thresh, boolean greaterThan){
      
      Img<T> out = img.factory().create(img, img.firstElement());
      RandomAccess<T> ra = out.randomAccess();
      Cursor<T> c = img.cursor();
      
      while(c.hasNext()){
         T t = c.next();
         ra.setPosition(c);
         
         if( greaterThan && t.getRealDouble() > thresh)
         {
            ra.get().set(t);
         }
         else if( !greaterThan && t.getRealDouble() < thresh )
         {
            ra.get().set(t);
         }
      }
      return out;
   }
   
   public static < T extends RealType< T >> Img<T> thresholdMap(Img<T> img, double thresh, boolean greaterThan){
	      
	      Img<T> out = img.factory().create(img, img.firstElement());
	      RandomAccess<T> ra = out.randomAccess();
	      Cursor<T> c = img.cursor();
	      
	      while(c.hasNext()){
	         T t = c.next();
	         ra.setPosition(c);
	         
	         if( greaterThan && t.getRealDouble() > thresh)
	         {
	            ra.get().setOne();
	         }
	         else if( !greaterThan && t.getRealDouble() < thresh )
	         {
	            ra.get().setOne();
	         }
	      }
	      return out;
	   }
   
   public static < T extends NativeType< T >> ImagePlusImg<T, ?> copyToImagePlus(Img<T> img) throws Exception {


      ImagePlusImgFactory<T> factory = new ImagePlusImgFactory<T>();
      T t = img.randomAccess().get();
      ImagePlusImg<T, ?> ipImg = factory.create(img,t);

      logger.debug("create image plus of type: " + t.getClass());
      logger.debug("result is of type: " + ipImg.firstElement().getClass());

      Cursor<T> c_in  = ipImg.cursor();
      RandomAccess<T> ra = img.randomAccess();

      while(c_in.hasNext()){
         c_in.fwd();
         ra.setPosition(c_in);
         ra.get().set(c_in.get());

      }
      
//      logger.debug("output ip has " + ImgUtil.numNonZero(ipImg) +  " non zero values");

      return ipImg;

   }
   
   public static < T extends NativeType< T >> void copyInto(
		   RandomAccessible<T> src, 
		   IterableInterval<T> dest)
   {
	   Cursor<T> c_in  = dest.cursor();
	   RandomAccess<T> ra = src.randomAccess();

	   while(c_in.hasNext()){
		   c_in.fwd();
		   ra.setPosition(c_in);
		   ra.get().set(c_in.get());
	   }
   }
   
   public static < T extends NativeType< T >> ImagePlusImg<T, ?> copyToImagePlus(
		   		RandomAccessible<T> img,
		   		IterableInterval<?> interval
		   ) throws Exception {


	      ImagePlusImgFactory<T> factory = new ImagePlusImgFactory<T>();
	      T t = img.randomAccess().get();
	      ImagePlusImg<T, ?> ipImg = factory.create(interval,t);

	      logger.debug("create image plus of type: " + t.getClass());
	      logger.debug("result is of type: " + ipImg.firstElement().getClass());

	      Cursor<T> c_in  = ipImg.cursor();
	      RandomAccess<T> ra = img.randomAccess();

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

   public static <L extends AbstractIntegerType<L>> ArrayList<Integer> uniqueInt( Img<L> img ){
	   
	   ArrayList<Integer> set = new ArrayList<Integer>();
	   Cursor<L> cursor = img.cursor();
	   while(cursor.hasNext()){
		   int l = cursor.next().getInteger();
		   if(!set.contains(l)){
			   set.add(l);
		   }
	   }
	   
	   return set;
   }
   
   public static HashSet<Float> unique( ImagePlus img ){
	   
	   HashSet<Float> set = new HashSet<Float>();
	   ImageProcessor ip = img.getProcessor();
	   int N  = ip.getPixelCount();
	   for (int i=0; i<N; i++) {
		   set.add(  ip.getf(i)  );
	   }
	   
	   return set;
   }
   
   public static <L extends AbstractIntegerType<L>> HashSet<L> unique2( Img<L> img ){
	   
	   HashSet<L> set = new HashSet<L>();
	   Cursor<L> cursor = img.cursor();
	   while(cursor.hasNext()){
		   set.add(cursor.next());
	   }
	   
	   return set;
   }
   
   public static <L extends AbstractIntegerType<L>> void combineValues( Img<L> img, int[][] spec, boolean strict){

      removeDuplicates(spec);

      Cursor<L> cursor = img.cursor();
      while(cursor.hasNext()){
         cursor.fwd();
         
         int currLab = cursor.get().getInteger(); // current label
         if(strict){
            cursor.get().setInteger( searchStrict( spec, currLab) );
         }else{
            cursor.get().setInteger( search( spec, currLab) );
         }

      }
   }

   public static <L extends AbstractIntegerType<L>> void replaceValues( Img<L> img, int[][] spec ){

      removeDuplicates(spec);

      Cursor<L> cursor = img.cursor();
      while(cursor.hasNext()){
         cursor.fwd();
         
         int currLab = cursor.get().getInteger(); // current label
         cursor.get().setInteger( replace( spec, currLab) );

      }
   }
  /** 
    * returns the (first) row in array that contains val
    * @param array
    * @param val 
    * @return
    */
   private static int searchStrict(int[][] array, int val){
      for(int i=0; i<array.length; i++){
         for(int j=0; j<array[0].length; j++){
            if(array[i][j]==val){
               return i;
            }   
         }   
      }   
      return -1; 
   }   

  /** 
    * returns the (first) row in array that contains val
    * @param array
    * @param val 
    * @return
    */
   private static int search(int[][] array, int val){
      for(int i=0; i<array.length; i++){
         for(int j=0; j<array[0].length; j++){
            if(array[i][j]==val){
               return i;
            }   
         }   
      }   
      return val; 
   }   
   
   /** 
    * returns the (first) row in array that contains val
    * @param array
    * @param val 
    * @return
    */
   private static int replace(int[][] array, int val){
      for(int i=0; i<array.length; i++){
         if(array[i][0]==val){
            return array[i][1];
         }   
      }   
      return val;
   }  

	private static void removeDuplicates(int[][] in){
		
		int currentVal = 0;
		for(int i = 0; i < in.length;i++ )
			for(int j = 0; j < in[i].length;j++){
				
				currentVal = in[i][j];
				if(currentVal > Integer.MIN_VALUE){
					for(int ii = i; ii < in.length;ii++ )
						for(int jj = 0; jj < in[ii].length;jj++){
							if(!(i == ii && j == jj) && in[ii][jj] == in[i][j]){
								in[ii][jj] = Integer.MIN_VALUE;
							}
						}
				}
			}
		
	}
	
	public <T extends RealType<T>> int numLessThanZero(Img<T> in) {
		int count = 0;
		Cursor<T> cursor = in.cursor();

		while (cursor.hasNext()) {
			if (cursor.next().getRealDouble() <= 0) {
				count++;
			}
		}
		return count;
	}

	public int numLessThanZero(float[][][] in) {
		int nx = in.length;
		int ny = in[0].length;
		int nz = in[0][0].length;

		int count = 0;

		for (int x = 0; x < nx; x++)
			for (int y = 0; y < ny; y++)
				for (int z = 0; z < nz; z++) {
					if (in[x][y][z] <= 0) {
						count++;
					}
				}
		return count;
	}

}
