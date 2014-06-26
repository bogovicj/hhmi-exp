package net.imglib2.util;

import net.imglib2.*;
import net.imglib2.img.*;
import net.imglib2.type.numeric.*;
import net.imglib2.view.Views;
import net.imglib2.algorithm.gauss3.*;
import net.imglib2.exception.IncompatibleTypeException;

public class Resampling {
	
	/**
	 * Create a downsampled {@link Img}.
	 *
	 * @param img the source image
	 * @param downsampleFactors scaling factors in each dimension
	 * @param sourceSigmas the Gaussian at which the source was sampled (guess 0.5 if you do not know)
	 * @param targetSigmas the Gaussian at which the target will be sampled
	 *
	 * @return a new {@link Img}
	 */
	public static <T extends RealType<T>> Img<T> resampleGaussian( 
			RandomAccessibleInterval<T> img,
			ImgFactory<T> factory,
			double[] downsampleFactors, 
			double[] sourceSigmas,
			double[] targetSigmas)
	{
		
		int ndims = img.numDimensions();
		long[] sz = new long[ ndims ];
		
		img.dimensions(sz);
		
		double[] sigs = new double[ ndims ];
		long[] newSz  = new long[ ndims ];
		
		for( int d = 0; d<ndims; d++)
		{
			double s = targetSigmas[d] * downsampleFactors[d]; 
			sigs[d] = Math.sqrt( s * s  - sourceSigmas[d] * sourceSigmas[d] );
			
			newSz[d] = Math.round( sz[d] / downsampleFactors[d] );
		}
		
		int[] pos = new int[ ndims ];
		
		RandomAccess<T> inRa = img.randomAccess();
		inRa.setPosition(pos);
		
		Img<T> out = factory.create( newSz, inRa.get());
		Img<T> tmp = factory.create( sz, inRa.get());
		
		try {
			Gauss3.gauss(sigs, Views.extendBorder(img), tmp);
		} catch (IncompatibleTypeException e) {
			e.printStackTrace();
		}
		
		
		int[][] coordLUT = new int[ndims][];
		for( int d = 0; d<ndims; d++)
		{
			int nd = (int)img.dimension(d);
			coordLUT[d] = new int[nd];
			
			for( int i = 0; i<nd; i++)
			{
				coordLUT[d][i] =  (int)Math.min( img.max(d), Math.max( 0, Math.round( i * downsampleFactors[d] ) ));
			}
		}

		Cursor<T> outc = out.cursor();
		
		while( outc.hasNext() )
		{
			outc.fwd();
			outc.localize( pos );
			setPositionFromLut( pos, coordLUT, inRa );
			outc.get().set( inRa.get() );
		}
		
		return out;
	}

	/**
	 * 
	 * @param destPos
	 * @param lut
	 * @param srcPos
	 */
	public static void setPositionFromLut( int[] destPos, int[][] lut, Positionable srcPos )
	{
		for ( int d = 0; d < destPos.length; d++ )
			srcPos.setPosition( lut[ d ][ destPos[d] ], d);
	}
	
}
