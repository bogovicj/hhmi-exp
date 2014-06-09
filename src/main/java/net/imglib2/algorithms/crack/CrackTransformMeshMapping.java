package net.imglib2.algorithms.crack;

import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.interpolation.InterpolatorFactory;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.view.Views;
import mpicbg.ij.InvertibleTransformMapping;
import mpicbg.models.AffineModel2D;
import mpicbg.models.PointMatch;
import mpicbg.util.Util;

public class CrackTransformMeshMapping< T extends CrackTransformMesh > extends InvertibleTransformMapping< T > 
{

	final static private class MapTriangleThread extends Thread
	{
		final private AtomicInteger i;
		final private List< AffineModel2D > triangles;
		final private CrackTransformMesh transform;
		final ImageProcessor source, target;
		MapTriangleThread(
				final AtomicInteger i,
				final List< AffineModel2D > triangles,
				final CrackTransformMesh transform,
				final ImageProcessor source,
				final ImageProcessor target )
		{
			this.i = i;
			this.triangles = triangles;
			this.transform = transform;
			this.source = source;
			this.target = target;
		}
		
		final public void run()
		{
			int k = i.getAndIncrement();
			while ( !isInterrupted() && k < triangles.size() )
			{
				mapTriangle( transform, triangles.get( k ), source, target );
				k = i.getAndIncrement();
			}
		}
	}
	
//	final static private class MapTriangleInterpolatedThread extends Thread
//	{
//		final private AtomicInteger i;
//		final private List< AffineModel2D > triangles;
//		final private CrackTransformMesh transform;
//		final ImageProcessor source, target;
//		MapTriangleInterpolatedThread(
//				final AtomicInteger i,
//				final List< AffineModel2D > triangles,
//				final CrackTransformMesh transform,
//				final ImageProcessor source,
//				final ImageProcessor target )
//		{
//			this.i = i;
//			this.triangles = triangles;
//			this.transform = transform;
//			this.source = source;
//			this.target = target;
//		}
//		
//		final public void run()
//		{
//			int k = i.getAndIncrement();
//			while ( !isInterrupted() && k < triangles.size() )
//			{
//				mapTriangleInterpolated( transform, triangles.get( k ), source, target );
//				k = i.getAndIncrement();
//			}
//		}
//	}
	
	final static private class MapTriangleInverseThread extends Thread
	{
		final private AtomicInteger i;
		final private List< AffineModel2D > triangles;
		final private CrackTransformMesh transform;
		final ImageProcessor source, target;
		MapTriangleInverseThread(
				final AtomicInteger i,
				final List< AffineModel2D > triangles,
				final CrackTransformMesh transform,
				final ImageProcessor source,
				final ImageProcessor target )
		{
			this.i = i;
			this.triangles = triangles;
			this.transform = transform;
			this.source = source;
			this.target = target;
		}
		
		final public void run()
		{
			int k = i.getAndIncrement();
			while ( !isInterrupted() && k < triangles.size() )
			{
				mapTriangleInverse( transform, triangles.get( k ), source, target );
				k = i.getAndIncrement();
			}
		}
	}
	
	final static private class MapTriangleInverseInterpolatedThread extends Thread
	{
		final private AtomicInteger i;
		final private List< AffineModel2D > triangles;
		final private CrackTransformMesh transform;
		final ImageProcessor source, target;
		MapTriangleInverseInterpolatedThread(
				final AtomicInteger i,
				final List< AffineModel2D > triangles,
				final CrackTransformMesh transform,
				final ImageProcessor source,
				final ImageProcessor target )
		{
			this.i = i;
			this.triangles = triangles;
			this.transform = transform;
			this.source = source;
			this.target = target;
		}
		
		final public void run()
		{
			int k = i.getAndIncrement();
			while ( !isInterrupted() && k < triangles.size() )
			{
				mapTriangleInverseInterpolated( transform, triangles.get( k ), source, target );
				k = i.getAndIncrement();
			}
		}
	}
	
	
	public CrackTransformMeshMapping( final T t )
	{
		super( t );
	}
	
	/**
	 * 
	 * @param pm PointMatches
	 * @param min x = min[0], y = min[1]
	 * @param max x = max[0], y = max[1]
	 */
	final static protected void calculateBoundingBox(
			final ArrayList< PointMatch > pm,
			final float[] min,
			final float[] max )
	{
		final float[] first = pm.get( 0 ).getP2().getW();
		min[ 0 ] = first[ 0 ];
		min[ 1 ] = first[ 1 ];
		max[ 0 ] = first[ 0 ];
		max[ 1 ] = first[ 1 ];
		
		for ( final PointMatch p : pm )
		{
			final float[] t = p.getP2().getW();
			if ( t[ 0 ] < min[ 0 ] ) min[ 0 ] = t[ 0 ];
			else if ( t[ 0 ] > max[ 0 ] ) max[ 0 ] = t[ 0 ];
			if ( t[ 1 ] < min[ 1 ] ) min[ 1 ] = t[ 1 ];
			else if ( t[ 1 ] > max[ 1 ] ) max[ 1 ] = t[ 1 ];
		}
	}
	
	
	/**
	 * Checks if a location is inside a given triangle.
	 * 
	 * @param pm
	 * @param t
	 * @return
	 */
	final static protected boolean isInTriangle(
			final float ax,
			final float ay,
			final float bx,
			final float by,
			final float cx,
			final float cy,
			final float tx,
			final float ty )
	{
		final boolean d;
		{
			final float x1 = bx - ax;
			final float y1 = by - ay;
			final float x2 = tx - ax;
			final float y2 = ty - ay;
			d = x1 * y2 - y1 * x2 < 0;
		}
		{
			final float x1 = cx - bx;
			final float y1 = cy - by;
			final float x2 = tx - bx;
			final float y2 = ty - by;
			if ( d ^ x1 * y2 - y1 * x2 < 0 ) return false;
		}
		{
			final float x1 = ax - cx;
			final float y1 = ay - cy;
			final float x2 = tx - cx;
			final float y2 = ty - cy;
			if ( d ^ x1 * y2 - y1 * x2 < 0 ) return false;
		}
		return true;
	}
	
	
	final static protected void mapTriangle(
			final CrackTransformMesh m, 
			final AffineModel2D ai,
			final ImageProcessor source,
			final ImageProcessor target )
	{
		final int w = target.getWidth() - 1;
		final int h = target.getHeight() - 1;
		final ArrayList< PointMatch > pm = m.getAV().get( ai );
		final float[] min = new float[ 2 ];
		final float[] max = new float[ 2 ];
		calculateBoundingBox( pm, min, max );
		final int minX = Math.max( 0, Util.roundPos( min[ 0 ] ) );
		final int minY = Math.max( 0, Util.roundPos( min[ 1 ] ) );
		final int maxX = Math.min( w, Util.roundPos( max[ 0 ] ) );
		final int maxY = Math.min( h, Util.roundPos( max[ 1 ] ) );
		
		final float[] a = pm.get( 0 ).getP2().getW();
		final float ax = a[ 0 ];
		final float ay = a[ 1 ];
		final float[] b = pm.get( 1 ).getP2().getW();
		final float bx = b[ 0 ];
		final float by = b[ 1 ];
		final float[] c = pm.get( 2 ).getP2().getW();
		final float cx = c[ 0 ];
		final float cy = c[ 1 ];
		final float[] t = new float[ 2 ];
		for ( int y = minY; y <= maxY; ++y )
		{
			for ( int x = minX; x <= maxX; ++x )
			{
				if ( isInTriangle( ax, ay, bx, by, cx, cy, x, y ) )
				{
					t[ 0 ] = x;
					t[ 1 ] = y;
					try
					{
						ai.applyInverseInPlace( t );
					}
					catch ( Exception e )
					{
						//e.printStackTrace( System.err );
						continue;
					}
					target.putPixel( x, y, source.getPixel( ( int )( t[ 0 ] + 0.5f ), ( int )( t[ 1 ] + 0.5f ) ) );
				}
			}
		}
	}
	
			
	final static protected void mapTriangleInterpolated(
			final CrackTransformMesh m, 
			final Collection<AffineModel2D> aiList,
			final ImageProcessor source,
			final ImageProcessor target )
	{
		final int w = target.getWidth() - 1;
		final int h = target.getHeight() - 1;

		final float[] t = new float[ 2 ];
		for ( int y = 0; y <= h; ++y )
		{
			for ( int x = 0; x <= w; ++x )
			{
				
				//boolean hitThis = false;
				boolean inTriangle = false;
				for( AffineModel2D ai : aiList ){
					
					ArrayList< PointMatch > pm = m.getAV().get( ai );
					
					float[] a = pm.get( 0 ).getP2().getW();
					float ax = a[ 0 ];
					float ay = a[ 1 ];
					float[] b = pm.get( 1 ).getP2().getW();
					float bx = b[ 0 ];
					float by = b[ 1 ];
					float[] c = pm.get( 2 ).getP2().getW();
					float cx = c[ 0 ];
					float cy = c[ 1 ];
					
					if ( isInTriangle( ax, ay, bx, by, cx, cy, x, y ) )
					{
						
						boolean isInCrack = m.getAC().get(ai).booleanValue();
//						if( isInCrack ){
//							System.out.println(" crack-tastic\n"  +
//							"("+ax + "," + ay + ")  (" + bx + ","+by+")  (" + cx+","+cy +")");
//						}
						
						t[ 0 ] = x;
						t[ 1 ] = y;
						try
						{
							ai.applyInverseInPlace( t );
						}
						catch ( Exception e )
						{
							//e.printStackTrace( System.err );
							System.out.println(" error at " + x + " " + y );
							continue;
						}
						
//						if( hitThis ){
//							System.out.println(" doubled up ");
//						}
						
						if( isInCrack ){
							target.putPixel( x, y, 0 );
							//hitThis = true;
						}else{
							target.putPixel( x, y, source.getPixelInterpolated( t[ 0 ], t[ 1 ] ) );
							//hitThis = true;
						}
						
						if( x == 109 && y ==95){
							ax += 0f;
						}
						
						inTriangle = true;
					}
				}
				
				if( !inTriangle ){
					target.putPixel( x, y, source.getPixel( x, y ) );
				}
			}
		}
	}
	
//	InterpolatorFactory<T,RandomAccessible<T>
	final static protected <T extends NumericType<T>> void mapTriangleInterpolated(
			final CrackTransformMesh m, 
			final Collection<AffineModel2D> aiList,
			final RandomAccessibleInterval<T> source,
			final IterableInterval<T> target )
	{
		NLinearInterpolatorFactory<T> interp = new NLinearInterpolatorFactory<T>(); 
		mapTriangleInterpolated(
				m, aiList, source,target, interp );
	}
	
	final static protected <T extends NumericType<T>> void mapTriangleInterpolated(
			final CrackTransformMesh m, 
			final Collection<AffineModel2D> aiList,
			final RandomAccessibleInterval<T> source,
			final IterableInterval<T> target,
			final InterpolatorFactory<T,RandomAccessible<T>> interp)
	{
		Cursor<T> curs = target.cursor();
		
		RealRandomAccessible<T> srcRa = Views.interpolate( 
					Views.extendZero( source ), 
					interp );
		RealRandomAccess<T> ra = srcRa.realRandomAccess();
		

		final float[] t = new float[ 2 ];
		while( curs.hasNext() )
		{
			curs.fwd();
			curs.localize( t );
			
			float x = t[0];
			float y = t[1];
			
			//boolean hitThis = false;
			boolean inTriangle = false;
			for( AffineModel2D ai : aiList ){

				ArrayList< PointMatch > pm = m.getAV().get( ai );

				float[] a = pm.get( 0 ).getP2().getW();
				float ax = a[ 0 ];
				float ay = a[ 1 ];
				float[] b = pm.get( 1 ).getP2().getW();
				float bx = b[ 0 ];
				float by = b[ 1 ];
				float[] c = pm.get( 2 ).getP2().getW();
				float cx = c[ 0 ];
				float cy = c[ 1 ];

				if ( isInTriangle( ax, ay, bx, by, cx, cy, x, y ) )
				{

					boolean isInCrack = m.getAC().get(ai).booleanValue();

					try
					{
						ai.applyInverseInPlace( t );
					}
					catch ( Exception e )
					{
						//e.printStackTrace( System.err );
						System.out.println(" error at " + x + " " + y );
						continue;
					}

					if( isInCrack ){
						curs.get().setZero();
						//hitThis = true;
					}else{
						ra.setPosition( t );
						curs.get().set( ra.get() );
						//hitThis = true;
					}

					inTriangle = true;
				}
			}

			if( !inTriangle ){
				try{
					ra.setPosition( t );
					curs.get().set( ra.get() );
				}catch( Exception e){
					System.out.println(" error at " + x + " " + y);
					continue;
				}
			}

		}
	}
	
	final static protected <T extends NumericType<T>> void mapMask(
			final CrackTransformMesh m, 
			final IterableInterval<T> target)
	{
		Cursor<T> curs = target.cursor();
		Set<AffineModel2D> aiList = m.av.keySet();

		final float[] t = new float[ 2 ];
		while( curs.hasNext() )
		{
			curs.fwd();
			curs.localize( t );
			
			float x = t[0];
			float y = t[1];
			
			for( AffineModel2D ai : aiList ){

				ArrayList< PointMatch > pm = m.getAV().get( ai );

				float[] a = pm.get( 0 ).getP2().getW();
				float ax = a[ 0 ];
				float ay = a[ 1 ];
				float[] b = pm.get( 1 ).getP2().getW();
				float bx = b[ 0 ];
				float by = b[ 1 ];
				float[] c = pm.get( 2 ).getP2().getW();
				float cx = c[ 0 ];
				float cy = c[ 1 ];

				if ( isInTriangle( ax, ay, bx, by, cx, cy, x, y ) )
				{

					boolean isInCrack = m.getAC().get(ai).booleanValue();

					try
					{
						ai.applyInverseInPlace( t );
					}
					catch ( Exception e )
					{
						//e.printStackTrace( System.err );
//						System.out.println(" error at " + x + " " + y );
						continue;
					}

					if( isInCrack ){
						curs.get().setOne();
						curs.get().mul( 255 );
					
					}
				}
			}
		}
	}
	
	
	final public void map(
			final ImageProcessor source,
			final ImageProcessor target,
			final int numThreads )
	{
		if ( numThreads == 1 )
		{
			/* no overhead for thread creation */
			final Set< AffineModel2D > s = transform.getAV().keySet();
			for ( final AffineModel2D ai : s )
				mapTriangle( transform, ai, source, target );
		}
		else
		{
			final List< AffineModel2D > l = new ArrayList< AffineModel2D >();
			l.addAll( transform.getAV().keySet() );
			final AtomicInteger i = new AtomicInteger( 0 );
			final ArrayList< Thread > threads = new ArrayList< Thread >( numThreads );
			for ( int k = 0; k < numThreads; ++k )
			{
				final Thread mtt = new MapTriangleThread( i, l, transform, source, target );
				threads.add( mtt );
				mtt.start();
			}
			for ( final Thread mtt : threads )
			{
				try
				{
					mtt.join();
				}
				catch ( InterruptedException e ) {}
			}
		}
	}
	
	@Override
	final public void map(
			final ImageProcessor source,
			final ImageProcessor target )
	{
		map( source, target, 1 );
	}
	
	final public void mapInterpolated(
			final ImageProcessor source,
			final ImageProcessor target,
			final int numThreads )
	{
		if ( numThreads == 1 )
		{
			/* no overhead for thread creation */
			final Set< AffineModel2D > s = transform.getAV().keySet();
			mapTriangleInterpolated( transform, s, source, target );
		}
		
	}
	
	@Override
	final public void mapInterpolated(
			final ImageProcessor source,
			final ImageProcessor target )
	{
		mapInterpolated( source, target, 1);
	}
	
	final public <F extends NumericType<F>> void mapInterpolated(
			final RandomAccessibleInterval<F> source,
			final IterableInterval<F> target )
	{
		final Set< AffineModel2D > s = transform.getAV().keySet();
		mapTriangleInterpolated( transform, s, source, target );
	}
	final public <F extends NumericType<F>> void mapInterpolated(
			final RandomAccessibleInterval<F> source,
			final IterableInterval<F> target,
			final InterpolatorFactory<F,RandomAccessible<F>> interp)
	{
		final Set< AffineModel2D > s = transform.getAV().keySet();
		mapTriangleInterpolated( transform, s, source, target, interp );
	}
	
	/**
	 * 
	 * @param pm PointMatches
	 * @param min x = min[0], y = min[1]
	 * @param max x = max[0], y = max[1]
	 */
	final static protected void calculateBoundingBoxInverse(
			final ArrayList< PointMatch > pm,
			final float[] min,
			final float[] max )
	{
		final float[] first = pm.get( 0 ).getP1().getL();
		min[ 0 ] = first[ 0 ];
		min[ 1 ] = first[ 1 ];
		max[ 0 ] = first[ 0 ];
		max[ 1 ] = first[ 1 ];
		
		for ( final PointMatch p : pm )
		{
			final float[] t = p.getP1().getL();
			if ( t[ 0 ] < min[ 0 ] ) min[ 0 ] = t[ 0 ];
			else if ( t[ 0 ] > max[ 0 ] ) max[ 0 ] = t[ 0 ];
			if ( t[ 1 ] < min[ 1 ] ) min[ 1 ] = t[ 1 ];
			else if ( t[ 1 ] > max[ 1 ] ) max[ 1 ] = t[ 1 ];
		}
	}
	
	final static protected void mapTriangleInverse(
			final CrackTransformMesh m, 
			final AffineModel2D ai,
			final ImageProcessor source,
			final ImageProcessor target )
	{
		final int w = target.getWidth() - 1;
		final int h = target.getHeight() - 1;
		final ArrayList< PointMatch > pm = m.getAV().get( ai );
		final float[] min = new float[ 2 ];
		final float[] max = new float[ 2 ];
		calculateBoundingBoxInverse( pm, min, max );
		final int minX = Math.max( 0, Util.roundPos( min[ 0 ] ) );
		final int minY = Math.max( 0, Util.roundPos( min[ 1 ] ) );
		final int maxX = Math.min( w, Util.roundPos( max[ 0 ] ) );
		final int maxY = Math.min( h, Util.roundPos( max[ 1 ] ) );
		
		final float[] a = pm.get( 0 ).getP1().getL();
		final float ax = a[ 0 ];
		final float ay = a[ 1 ];
		final float[] b = pm.get( 1 ).getP1().getL();
		final float bx = b[ 0 ];
		final float by = b[ 1 ];
		final float[] c = pm.get( 2 ).getP1().getL();
		final float cx = c[ 0 ];
		final float cy = c[ 1 ];
		final float[] t = new float[ 2 ];
		for ( int y = minY; y <= maxY; ++y )
		{
			for ( int x = minX; x <= maxX; ++x )
			{
				if ( isInTriangle( ax, ay, bx, by, cx, cy, x, y ) )
				{
					t[ 0 ] = x;
					t[ 1 ] = y;
					ai.applyInPlace( t );
					target.putPixel( x, y, source.getPixel( ( int )( t[ 0 ] + 0.5f ), ( int )( t[ 1 ] + 0.5f ) ) );
				}
			}
		}
	}
	
	final static protected void mapTriangleInverseInterpolated(
			final CrackTransformMesh m, 
			final AffineModel2D ai,
			final ImageProcessor source,
			final ImageProcessor target )
	{
		final int w = target.getWidth() - 1;
		final int h = target.getHeight() - 1;
		final ArrayList< PointMatch > pm = m.getAV().get( ai );
		final float[] min = new float[ 2 ];
		final float[] max = new float[ 2 ];
		calculateBoundingBoxInverse( pm, min, max );
		final int minX = Math.max( 0, Util.roundPos( min[ 0 ] ) );
		final int minY = Math.max( 0, Util.roundPos( min[ 1 ] ) );
		final int maxX = Math.min( w, Util.roundPos( max[ 0 ] ) );
		final int maxY = Math.min( h, Util.roundPos( max[ 1 ] ) );
		
		final float[] a = pm.get( 0 ).getP1().getL();
		final float ax = a[ 0 ];
		final float ay = a[ 1 ];
		final float[] b = pm.get( 1 ).getP1().getL();
		final float bx = b[ 0 ];
		final float by = b[ 1 ];
		final float[] c = pm.get( 2 ).getP1().getL();
		final float cx = c[ 0 ];
		final float cy = c[ 1 ];
		final float[] t = new float[ 2 ];
		for ( int y = minY; y <= maxY; ++y )
		{
			for ( int x = minX; x <= maxX; ++x )
			{
				if ( isInTriangle( ax, ay, bx, by, cx, cy, x, y ) )
				{
					t[ 0 ] = x;
					t[ 1 ] = y;
					ai.applyInPlace( t );
					target.putPixel( x, y, source.getPixelInterpolated( t[ 0 ], t[ 1 ] ) );
				}
			}
		}
	}
	
	final public void mapInverse(
			final ImageProcessor source,
			final ImageProcessor target,
			final int numThreads )
	{
		if ( numThreads == 1 )
		{
			/* no overhead for thread creation */
			final Set< AffineModel2D > s = transform.getAV().keySet();
			for ( final AffineModel2D ai : s )
				mapTriangleInverse( transform, ai, source, target );
		}
		else
		{
			final List< AffineModel2D > l = new ArrayList< AffineModel2D >();
			l.addAll( transform.getAV().keySet() );
			final AtomicInteger i = new AtomicInteger( 0 );
			final ArrayList< Thread > threads = new ArrayList< Thread >( numThreads );
			for ( int k = 0; k < numThreads; ++k )
			{
				final Thread mtt = new MapTriangleInverseThread( i, l, transform, source, target );
				threads.add( mtt );
				mtt.start();
			}
			for ( final Thread mtt : threads )
			{
				try
				{
					mtt.join();
				}
				catch ( InterruptedException e ) {}
			}
		}
	}
	
	@Override
	final public void mapInverse(
			final ImageProcessor source,
			final ImageProcessor target )
	{
		mapInverse( source, target, Runtime.getRuntime().availableProcessors() );
	}
	
	final public void mapInverseInterpolated(
			final ImageProcessor source,
			final ImageProcessor target,
			final int numThreads )
	{
		if ( numThreads == 1 )
		{
			/* no overhead for thread creation */
			final Set< AffineModel2D > s = transform.getAV().keySet();
			for ( final AffineModel2D ai : s )
				mapTriangleInverseInterpolated( transform, ai, source, target );
		}
		else
		{
			final List< AffineModel2D > l = new ArrayList< AffineModel2D >();
			l.addAll( transform.getAV().keySet() );
			final AtomicInteger i = new AtomicInteger( 0 );
			final ArrayList< Thread > threads = new ArrayList< Thread >( numThreads );
			for ( int k = 0; k < numThreads; ++k )
			{
				final Thread mtt = new MapTriangleInverseInterpolatedThread( i, l, transform, source, target );
				threads.add( mtt );
				mtt.start();
			}
			for ( final Thread mtt : threads )
			{
				try
				{
					mtt.join();
				}
				catch ( InterruptedException e ) {}
			}
		}
	}
	
	@Override
	final public void mapInverseInterpolated(
			final ImageProcessor source,
			final ImageProcessor target )
	{
		mapInverseInterpolated( source, target, Runtime.getRuntime().availableProcessors() );
	}
	
}

