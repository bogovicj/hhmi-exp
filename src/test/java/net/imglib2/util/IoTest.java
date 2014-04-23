package net.imglib2.util;

import static org.junit.Assert.*;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.integer.ByteType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.integer.ShortType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;

import org.junit.Test;

public class IoTest {

	@Test
	public void testConvertToImagPlus(){


		int[] sz = new int[]{32,32,32};

		ArrayImgFactory<ByteType> factory = new ArrayImgFactory<ByteType>();
		Img<ByteType> ubimg = factory.create( sz, new ByteType()); 
		
		RandomAccess<ByteType> ra = ubimg.randomAccess();
		ra.setPosition(new int[]{1,1,0});
		ra.get().set((byte) 4);
		
//		ImgUtil.write(ubimg, "/groups/jain/home/bogovicj/test/imglib2Io/ubyteImg.tif");
		
		
		Img<FloatType> gradimgF = ImgUtil.createGradientImgX(sz[0], sz[1], sz[2], new FloatType());
		ImgUtil.writeFloat(gradimgF, "/groups/jain/home/bogovicj/test/imglib2Io/floatGradX.tif");
		
		
		Img<UnsignedByteType> gradimgB = ImgUtil.createGradientImgX(sz[0], sz[1], sz[2], new UnsignedByteType());
		ImgUtil.writeByte(gradimgB, "/groups/jain/home/bogovicj/test/imglib2Io/byteGradX.tif");
		
		Img<ShortType> gradimgS = ImgUtil.createGradientImgX(sz[0], sz[1], sz[2], new ShortType());
		ImgUtil.writeShort(gradimgS, "/groups/jain/home/bogovicj/test/imglib2Io/shortGradX.tif");
		
	}

}
