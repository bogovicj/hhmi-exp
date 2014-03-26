package net.imglib2.util;

import edu.jhu.ece.iacl.utility.ArrayUtil;
import ij.ImagePlus;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.imageplus.ImagePlusImg;
import net.imglib2.img.imageplus.ImagePlusImgFactory;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.RealType;

public class ImgUtil {

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
}
