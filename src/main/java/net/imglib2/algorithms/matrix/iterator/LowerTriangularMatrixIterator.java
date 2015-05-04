package net.imglib2.algorithms.matrix.iterator;

import java.util.Iterator;

/**
 * 
 * @author John Bogovic <bogovicj@janelia.hhmi.org>
 * 
 */
public class LowerTriangularMatrixIterator implements Iterator< Long >
{

	protected long nrows; // width and height of matrix
	protected long numElements;

	protected long start; // starting point
	protected long n; // current index
	protected long count; // count
	protected long i, j; // current row and column in matrix

	public LowerTriangularMatrixIterator( int nrowsIn, int start )
	{
		this.nrows = nrowsIn;

		if ( nrows % 2 == 0 )
		{
			numElements = (nrows - 1) * (nrows / 2);
		} else
		{
			numElements = nrows * ((nrows - 1) / 2);
		}

		setStart( start );
	}

	public long getNum()
	{
		return numElements;
	}

	public LowerTriangularMatrixIterator( int nrows )
	{
		this( nrows, 1 ); // the "one" index is the
	}

	protected void setStart( int start )
	{
		this.start = start;
		n = start;
		long[] tmp = idxToCoords( start );
		i = tmp[ 0 ];
		j = tmp[ 1 ];
	}

	@Override
	public boolean hasNext()
	{
		return (count < numElements);
	}

	public long[] next( int K )
	{
		long[] nextN = new long[ K ];
		for ( int i = 0; i < K; i++ )
		{
			nextN[ i ] = next();
		}
		return nextN;
	}

	@Override
	public Long next()
	{
		long nOut = n; // the value to return
		count++;

		// for column major
		if ( i == nrows - 1 )
		{
			j++; // jump to the next column
			i = j + 1; // hop to the first column after the diagonal
			n += j + 2; // update the linear index

		} else if ( i < nrows )
		{
			i++;
			n++;
		}

		return nOut;
	}

	public long[] nextCoord()
	{
		return idxToCoords( next() );
	}

	public long[][] nextCoords( int K )
	{
		long[][] nextNcoords = new long[ K ][ 2 ];
		for ( int i = 0; i < K; i++ )
		{
			nextNcoords[ i ] = nextCoord();
		}
		return nextNcoords;
	}

	public void offset( int offset )
	{
		for ( int i = 0; i < offset; i++ )
		{
			if ( hasNext() )
				next();
		}

	}

	@Override
	public void remove()
	{
		// nothing to do
	}

	/**
	 * Column major order
	 * 
	 * @param i
	 *            row
	 * @param j
	 *            column
	 * @return the linear index
	 */
	public long coordsToIdx( int i, int j )
	{
		return i + nrows * j;
	}

	/**
	 * 
	 * @param n
	 *            the linear index
	 * @return array containing the { rowIndex, columnIndex }
	 */
	public long[] idxToCoords( long n )
	{
		long[] coords = new long[ 2 ];
		idxToCoords( n, coords );
		return coords;
	}

	/**
	 * 
	 * @param n
	 *            the linear indices
	 * @return array containing the { rowIndex, columnIndex }
	 */
	public long[][] idxToCoords( long[] n )
	{
		int K = n.length;
		long[][] coords = new long[ K ][ 2 ];
		for ( int i = 0; i < K; i++ )
		{
			idxToCoords( n[ i ], coords[ i ] );
		}
		return coords;
	}

	/**
	 * 
	 * @param n
	 *            the linear index
	 * @return array containing the { rowIndex, columnIndex }
	 */
	public void idxToCoords( long n, long[] coords )
	{

		coords[ 0 ] = (n % nrows);

		n -= coords[ 0 ];
		n /= nrows;

		coords[ 1 ] = (n % nrows);
	}

	public static void main( String[] args )
	{
		LowerTriangularMatrixIterator ltmi = new LowerTriangularMatrixIterator(
				50000 );

		System.out.println( "numel: " + ltmi.numElements );
		System.out.println( "num:" + ltmi.getNum() );

		int sm = 50000;
		int bg = 2 * sm;
		LowerTriangularMatrixIterator miSm = new LowerTriangularMatrixIterator(
				sm );
		LowerTriangularMatrixIterator miBg = new LowerTriangularMatrixIterator(
				bg );

		System.out.println( "num sm:" + miSm.getNum() );
		System.out.println( "num bg:" + miBg.getNum() );
		System.out.println( "is numSm < numBg:"
				+ (miSm.getNum() < miBg.getNum()) );

		// System.out.println( "(0,1): " + utmi.coordsToIdx( 0, 1 ) );
		// System.out.println( "(1,0): " + utmi.coordsToIdx( 1, 0 ) );
		//
		// System.out.println( "(0): "
		// + ArrayUtil.printArray( utmi.idxToCoords( 0 ) ) );
		// System.out.println( "(1): "
		// + ArrayUtil.printArray( utmi.idxToCoords( 1 ) ) );
		// System.out.println( "(5): "
		// + ArrayUtil.printArray( utmi.idxToCoords( 5 ) ) );

		/*
		 * while ( ltmi.hasNext() ) { System.out.println( ltmi.next() ); }
		 */

		System.out.println( "done" );
	}
}
