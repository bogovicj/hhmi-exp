/*
 * #%L
 * ImgLib2: a general-purpose, multidimensional image processing library.
 * %%
 * Copyright (C) 2009 - 2014 Stephan Preibisch, Tobias Pietzsch, Barry DeZonia,
 * Stephan Saalfeld, Albert Cardona, Curtis Rueden, Christian Dietz, Jean-Yves
 * Tinevez, Johannes Schindelin, Lee Kamentsky, Larry Lindsey, Grant Harris,
 * Mark Hiner, Aivar Grislis, Martin Horn, Nick Perry, Michael Zinsmaier,
 * Steffen Jaensch, Jan Funke, Mark Longair, and Dimiter Prodanov.
 * %%
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * #L%
 */

package net.imglib2.algorithms.region.localneighborhood;

import java.util.Iterator;

import net.imglib2.AbstractEuclideanSpace;
import net.imglib2.AbstractLocalizable;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.IterableRealInterval;
import net.imglib2.Positionable;
import net.imglib2.RandomAccess;
import net.imglib2.RealPositionable;

public class CrossNeighborhoodSkipCenter< T > extends AbstractLocalizable implements Neighborhood< T >
{
	public static < T > CrossNeighborhoodFactory< T > factory()
	{
		return new CrossNeighborhoodFactory< T >() {
			@Override
			public Neighborhood< T > create( final long[] position, final long[] currentMin, final long[] currentMax, final Interval span, final RandomAccess< T > sourceRandomAccess )
			{
				return new CrossNeighborhoodSkipCenter< T >( position, currentMin, currentMax, span, sourceRandomAccess );
			}
		};
	}

	protected final long[] currentMin;

	protected final long[] currentMax;

	protected final long[] dimensions;

	protected final RandomAccess< T > sourceRandomAccess;

	protected final Interval structuringElementBoundingBox;

	protected long size;

	CrossNeighborhoodSkipCenter( final long[] position, final long[] currentMin, final long[] currentMax, final Interval span, final RandomAccess< T > sourceRandomAccess )
	{
		super( position );
		this.currentMin = currentMin;
		this.currentMax = currentMax;
		dimensions = new long[ n ];
		span.dimensions( dimensions );

		size = dimensions[ 0 ];
		for ( int d = 1; d < n; ++d )
			size *= dimensions[ d ];

		this.sourceRandomAccess = sourceRandomAccess;
		this.structuringElementBoundingBox = span;
	}

	@Override
	public Interval getStructuringElementBoundingBox()
	{
		return structuringElementBoundingBox;
	}

	@Override
	public long size()
	{
		return size; // -1 because we skip the center pixel
	}

	@Override
	public T firstElement()
	{
		return cursor().next();
	}

	@Override
	public Object iterationOrder()
	{
		return this; // iteration order is only compatible with ourselves
	}

	@Override
	public boolean equalIterationOrder( final IterableRealInterval< ? > f )
	{
		return iterationOrder().equals( f.iterationOrder() );
	}

	@Override
	public double realMin( final int d )
	{
		return currentMin[ d ];
	}

	@Override
	public void realMin( final double[] min )
	{
		for ( int d = 0; d < n; ++d )
			min[ d ] = currentMin[ d ];
	}

	@Override
	public void realMin( final RealPositionable min )
	{
		for ( int d = 0; d < n; ++d )
			min.setPosition( currentMin[ d ], d );
	}

	@Override
	public double realMax( final int d )
	{
		return currentMax[ d ];
	}

	@Override
	public void realMax( final double[] max )
	{
		for ( int d = 0; d < n; ++d )
			max[ d ] = currentMax[ d ];
	}

	@Override
	public void realMax( final RealPositionable max )
	{
		for ( int d = 0; d < n; ++d )
			max.setPosition( currentMax[ d ], d );
	}

	@Override
	public Iterator< T > iterator()
	{
		return cursor();
	}

	@Override
	public long min( final int d )
	{
		return currentMin[ d ];
	}

	@Override
	public void min( final long[] min )
	{
		for ( int d = 0; d < n; ++d )
			min[ d ] = currentMin[ d ];
	}

	@Override
	public void min( final Positionable min )
	{
		for ( int d = 0; d < n; ++d )
			min.setPosition( currentMin[ d ], d );
	}

	@Override
	public long max( final int d )
	{
		return currentMax[ d ];
	}

	@Override
	public void max( final long[] max )
	{
		for ( int d = 0; d < n; ++d )
			max[ d ] = currentMax[ d ];
	}

	@Override
	public void max( final Positionable max )
	{
		for ( int d = 0; d < n; ++d )
			max.setPosition( currentMax[ d ], d );
	}

	@Override
	public void dimensions( final long[] dimensions )
	{
		for ( int d = 0; d < n; ++d )
			dimensions[ d ] = this.dimensions[ d ];
	}

	@Override
	public long dimension( final int d )
	{
		return dimensions[ d ];
	}

	@Override
	public LocalCrossCursor cursor()
	{
		return new LocalCrossCursor( sourceRandomAccess.copyRandomAccess() );
	}

	@Override
	public LocalCrossCursor localizingCursor()
	{
		return cursor();
	}

	public class LocalCrossCursor extends AbstractEuclideanSpace implements Cursor< T >
	{
		protected final RandomAccess< T > source;

		protected long index;
		
		protected int curdim;

		protected long maxIndexOnLine;

		public LocalCrossCursor( final RandomAccess< T > source )
		{
			super( source.numDimensions() );
			this.source = source;
			reset();
		}

		protected LocalCrossCursor( final LocalCrossCursor c )
		{
			super( c.numDimensions() );
			source = c.source.copyRandomAccess();
			index = c.index;
		}

		@Override
		public T get()
		{
			return source.get();
		}

		@Override
		public void fwd()
		{
			if ( ++index > maxIndexOnLine ){
				curdim++;
				nextLine();
				index = 0;
			}else{
				source.fwd( curdim );
			}
			
			// are we at the midpoint?
			boolean ismid = true;
			for ( int d = 0; d < n; ++d ){
				if(source.getIntPosition(d) != position[d]){
					ismid = false;
					break;
				}
			}
			
			// skip the midpoint
			if(ismid){ fwd(); }
			
		}

		private void nextLine()
		{
			for ( int d = 0; d < n; ++d )
			{
				if(d==curdim){
					source.setPosition(currentMin[d], d);
				}else{
					source.setPosition(position[d], d);
				}
			}
		}

		@Override
		public void jumpFwd( final long steps )
		{
			for ( long i = 0; i < steps; ++i )
				fwd();
		}

		@Override
		public T next()
		{
			fwd();
			return get();
		}

		@Override
		public void remove()
		{
			// NB: no action.
		}

		@Override
		public void reset()
		{
			index = 0;
			maxIndexOnLine = dimensions[ 0 ];
			curdim = 0;
			nextLine();
			source.bck( 0 );
		}

		@Override
		public boolean hasNext()
		{
			return (curdim < dimensions.length - 1) || (index < maxIndexOnLine - 1); 
		}

		@Override
		public float getFloatPosition( final int d )
		{
			return source.getFloatPosition( d );
		}

		@Override
		public double getDoublePosition( final int d )
		{
			return source.getDoublePosition( d );
		}

		@Override
		public int getIntPosition( final int d )
		{
			return source.getIntPosition( d );
		}

		@Override
		public long getLongPosition( final int d )
		{
			return source.getLongPosition( d );
		}

		@Override
		public void localize( final long[] position )
		{
			source.localize( position );
		}

		@Override
		public void localize( final float[] position )
		{
			source.localize( position );
		}

		@Override
		public void localize( final double[] position )
		{
			source.localize( position );
		}

		@Override
		public void localize( final int[] position )
		{
			source.localize( position );
		}

		@Override
		public LocalCrossCursor copy()
		{
			return new LocalCrossCursor( this );
		}

		@Override
		public LocalCrossCursor copyCursor()
		{
			return copy();
		}
	}
}
