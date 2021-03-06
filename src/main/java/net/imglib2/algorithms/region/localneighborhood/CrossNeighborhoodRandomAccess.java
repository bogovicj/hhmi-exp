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

import net.imglib2.Interval;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.region.localneighborhood.Neighborhood;

public final class CrossNeighborhoodRandomAccess< T > extends CrossNeighborhoodLocalizableSampler< T > implements RandomAccess< Neighborhood< T > >
{
	public CrossNeighborhoodRandomAccess( final RandomAccessibleInterval< T > source, final Interval span, final CrossNeighborhoodFactory< T > factory )
	{
		super( source, span, factory );
	}

	private CrossNeighborhoodRandomAccess( final CrossNeighborhoodRandomAccess< T > c )
	{
		super( c );
	}

	
	public void fwd( final int d )
	{
		++currentPos[ d ];
		++currentMin[ d ];
		++currentMax[ d ];
	}

	
	public void bck( final int d )
	{
		--currentPos[ d ];
		--currentMin[ d ];
		--currentMax[ d ];
	}

	
	public void move( final int distance, final int d )
	{
		currentPos[ d ] += distance;
		currentMin[ d ] += distance;
		currentMax[ d ] += distance;
	}

	
	public void move( final long distance, final int d )
	{
		currentPos[ d ] += distance;
		currentMin[ d ] += distance;
		currentMax[ d ] += distance;
	}

	
	public void move( final Localizable localizable )
	{
		for ( int d = 0; d < n; ++d )
		{
			final long distance = localizable.getLongPosition( d );
			currentPos[ d ] += distance;
			currentMin[ d ] += distance;
			currentMax[ d ] += distance;
		}
	}

	
	public void move( final int[] distance )
	{
		for ( int d = 0; d < n; ++d )
		{
			currentPos[ d ] += distance[ d ];
			currentMin[ d ] += distance[ d ];
			currentMax[ d ] += distance[ d ];
		}
	}

	
	public void move( final long[] distance )
	{
		for ( int d = 0; d < n; ++d )
		{
			currentPos[ d ] += distance[ d ];
			currentMin[ d ] += distance[ d ];
			currentMax[ d ] += distance[ d ];
		}
	}

	
	public void setPosition( final Localizable localizable )
	{
		for ( int d = 0; d < n; ++d )
		{
			final long position = localizable.getLongPosition( d );
			currentPos[ d ] = position;
			currentMin[ d ] = position + span.min( d );
			currentMax[ d ] = position + span.max( d );
		}
	}

	
	public void setPosition( final int[] position )
	{
		for ( int d = 0; d < n; ++d )
		{
			currentPos[ d ] = position[ d ];
			currentMin[ d ] = position[ d ] + span.min( d );
			currentMax[ d ] = position[ d ] + span.max( d );
		}
	}

	
	public void setPosition( final long[] position )
	{
		for ( int d = 0; d < n; ++d )
		{
			currentPos[ d ] = position[ d ];
			currentMin[ d ] = position[ d ] + span.min( d );
			currentMax[ d ] = position[ d ] + span.max( d );
		}
	}

	
	public void setPosition( final int position, final int d )
	{
		currentPos[ d ] = position;
		currentMin[ d ] = position + span.min( d );
		currentMax[ d ] = position + span.max( d );
	}

	
	public void setPosition( final long position, final int d )
	{
		currentPos[ d ] = position;
		currentMin[ d ] = position + span.min( d );
		currentMax[ d ] = position + span.max( d );
	}

	
	public CrossNeighborhoodRandomAccess< T > copy()
	{
		return new CrossNeighborhoodRandomAccess< T >( this );
	}

	
	public CrossNeighborhoodRandomAccess< T > copyRandomAccess()
	{
		return copy();
	}
}
