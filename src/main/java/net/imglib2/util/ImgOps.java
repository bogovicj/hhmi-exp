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
import net.imglib2.algorithms.patch.PatchTools;
import net.imglib2.exception.ImgLibException;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.imageplus.ImagePlusImg;
import net.imglib2.img.imageplus.ImagePlusImgFactory;
import net.imglib2.iterator.IntervalIterator;
import net.imglib2.meta.IntervalUtils;
import net.imglib2.ops.operation.randomaccessibleinterval.unary.DistanceMap;
import net.imglib2.type.NativeType;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.AbstractIntegerType;
import net.imglib2.type.numeric.integer.GenericByteType;
import net.imglib2.type.numeric.integer.GenericShortType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
import net.imglib2.view.composite.CompositeIntervalView;
import net.imglib2.view.composite.RealComposite;

public class ImgOps {

	public static Logger logger = LogManager.getLogger(ImgOps.class.getName());
	
	public static <T extends NativeType<T>> void write(Img<T> img, String fn)
	{
		try
		{
			ImagePlus ipdp = ImgOps.copyToImagePlus(img).getImagePlus();
			IJ.save(ipdp, fn);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	
	public static <T  extends RealType<T> & NativeType<T>> void writeFloat(Img<T> img, String fn)
	{
		logger.debug(" write float ");
		if( img instanceof ImagePlusImg ){
			logger.debug(" ij save ");
			try {
				IJ.save(((ImagePlusImg) img).getImagePlus(), fn);
			} catch (ImgLibException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			logger.debug(" save done");
			return;
		}
		try
		{
			ImagePlus ipdp = mimicImagePlusFloat( img, "name" );
			
			if(img.numDimensions()==2){
				copyToImageProcessor2dFloat(img, ipdp.getProcessor());
			}else if(img.numDimensions()==3){
				copyToImageProcessor3dFloat(img, ipdp);
			}
		
			logger.debug(" ij save ");
			IJ.save(ipdp, fn);
			logger.debug(" save done");
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
	public static <T extends RealType<T> & NativeType<T>> ImagePlus mimicImagePlusFloat(Img<T> img, String name)
	{
		logger.debug(" mimic image plus float ");
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
		logger.debug(" done mimic ");
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
	public static <T extends RealType<T> & NativeType<T>> void copyToImageProcessor2dFloat(Img<T> img, ImageProcessor ip) throws Exception {
		
		Cursor<T> c_in  = img.cursor();
		
		while(c_in.hasNext()){
			c_in.fwd();
			
			ip.setf(c_in.getIntPosition(0), c_in.getIntPosition(1), 
					c_in.get().getRealFloat());
			
		}
	}
	
	public static <T extends RealType<T> & NativeType<T>> void copyToImageProcessor3dFloat(Img<T> img, ImagePlus ip) throws Exception {
		
		Cursor<T> c_in  = img.cursor();
		logger.debug("copying to float 3d");
		int i = 0;
		while(c_in.hasNext()){
			c_in.fwd();
			
//			if ( i%1000 == 0 ) {logger.debug(" c_in " + (i++) ); }
//			else{ i++; }
			
			ImageProcessor fp = ip.getStack().getProcessor(c_in.getIntPosition(2)+1);
			
			fp.setf(c_in.getIntPosition(0), c_in.getIntPosition(1), 
					c_in.get().getRealFloat());
			
		}
		logger.debug("done copying");
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
	
	public static <T extends RealType<T>> Img<T> collapseSum(Img<T> in){
		
		RandomAccess<T> inra = in.randomAccess();
		
		long[] indim = IntervalUtils.getDims(in);
		long[] outdim = new long[ indim.length - 1 ];
		for (int i=0; i<outdim.length; i++){
			outdim[i] = indim[i];
		}
		int N = in.numDimensions();
		
		Img<T> out = in.factory().create(outdim, in.firstElement());
		RandomAccess<T> outra = out.randomAccess();
		
		IntervalIterator c = new IntervalIterator( outdim );
		while( c.hasNext() )
		{
			c.fwd();
			outra.setPosition(c);
			
			for (int i=0; i<outdim.length; i++){
				inra.setPosition( c.getIntPosition(i), i);
			}
			
			for ( int n=0; n < in.dimension(N-1); n++ )
			{
				inra.setPosition(n, (N-1));
				outra.get().add( inra.get() );
			}
		}
		
		return out;
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
		Img<B> maskInv = ImgOps.thresholdMap(mask, 1, false);
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
   
   public static <T extends NativeType<T> & RealType<T>> Img<T> createCheckerImg(int[] sz,  T t, int numLevels ){

	   if( numLevels < 2 )
	   {
		   logger.error("a checkered image can only be created with two or more levels - returning null");
		   return null;
	   }
	   
	   ArrayImgFactory<T> factory = new ArrayImgFactory<T>();
	   Img<T> out = factory.create( sz, t);
	   
	   Cursor<T> c = out.cursor();
	   int[] pos = new int[sz.length];
	   while ( c.hasNext() )
	   {
		   c.fwd();
		   c.localize(pos);
		   
		   int sum = ArrayUtil.sum(pos);
		   c.get().setReal(
				   	( sum % numLevels)
				   );
	   }
	   
	   return out;
   }

   public static <T extends NativeType<T> & RealType<T>> Img<T> createEdgeImg(int[] sz, double[] w,  T t, double sigma){
      ArrayImgFactory<T> factory = new ArrayImgFactory<T>();
      Img<T> out = factory.create( sz, t);
      
      double[] ctr = ArrayUtil.toDouble(sz);
      double[] ones = new double[ctr.length]; 
      ArrayUtil.fill(ones, 1);
      ctr = ArrayUtil.subtract(ctr, ones);
      ArrayUtil.divide(ctr, 2);
      
      logger.debug(" ctr = " + ArrayUtil.printArray(ctr));
      
      Cursor<T> c = out.localizingCursor();
      double[] pos = new double[out.numDimensions()];
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
   
   public static <T extends NativeType<T> & RealType<T>> Img<T> createGradientImgX( int[] sz, T t){
      ArrayImgFactory<T> factory = new ArrayImgFactory<T>();
      Img<T> out = factory.create(sz, t);
      
      Cursor<T> c = out.localizingCursor();
      int[] pos = new int[3];
      while(c.hasNext()){
         T val = c.next();
         c.localize(pos);
         val.setReal(pos[0]);
      }
      
      return out;
   }

   public static <T extends NativeType<T> & RealType<T>> Img<T> createGradientImgY( int[] sz, T t){
      ArrayImgFactory<T> factory = new ArrayImgFactory<T>();
      Img<T> out = factory.create(sz, t);
      
      Cursor<T> c = out.localizingCursor();
      int[] pos = new int[3];
      while(c.hasNext()){
         T val = c.next();
         c.localize(pos);
         val.setReal(pos[1]);
      }
      
      return out;
   }
   
   public static <T extends NativeType<T> & RealType<T>> Img<T> createGradientImgZ( int[] sz, T t){
      ArrayImgFactory<T> factory = new ArrayImgFactory<T>();
      Img<T> out = factory.create( sz, t);
      
      Cursor<T> c = out.localizingCursor();
      int[] pos = new int[3];
      while(c.hasNext()){
         T val = c.next();
         c.localize(pos);
         val.setReal(pos[2]);
      }
      
      return out;
   }
   
   public static <T extends NativeType<T> & RealType<T>> Img<T> createGaussianEllipseImg(int[] size, double[] ctr, double[] sigmas, double min, double max, T t){
	   return createGaussianEllipseImg( new ArrayImgFactory<T>(),
			   size, ctr, sigmas, min, max, t);
   }
   
   public static <T extends NativeType<T> & RealType<T>> Img<T> createGaussianEllipseImg(ImgFactory<T> factory, int[] size, double[] ctr, double[] sigmas, double min, double max, T t){
	      
	      Img<T> out = factory.create(size, t);
	      
	      Cursor<T> c = out.localizingCursor();
	      double[] pos = new double[3];
	      
	      while(c.hasNext())
	      {
	    	 c.fwd();
	         c.localize(pos);
	         
	         double res = 0;
	         for(int d=0; d<size.length; d++)
	         {
	        	 res +=  (pos[d] - ctr[d]) * (pos[d] - ctr[d]) / sigmas[d];
	         }
	         res = 1/(1 + res);
	         
	         if( res > min && res < max)  c.get().setReal(res);
	         
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
   
   public static < A extends RealType<A>, B extends RealType<B>> void copyInto(
		   RandomAccessible<A> src, 
		   IterableInterval<B> dest)
   {
	   Cursor<B> c_dest  = dest.cursor();
	   RandomAccess<A> ra = src.randomAccess();

	   while(c_dest.hasNext()){
		   c_dest.fwd();
		   ra.setPosition(c_dest);
		   c_dest.get().setReal( ra.get().getRealDouble() );
//		   System.out.println(" ra.get: " + ra.get());
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
	
	public static <T> void printPosition(RandomAccess<T> ra){
		int N = ra.numDimensions();
		for(int d=0; d<N; d++){
			System.out.print( " " + ra.getLongPosition(d));
		}
		System.out.print("\n");
	}

}
