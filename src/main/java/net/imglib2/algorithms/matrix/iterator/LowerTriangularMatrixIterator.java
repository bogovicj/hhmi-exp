package net.imglib2.algorithms.matrix.iterator;

import java.util.Iterator;

/**
 * 
 * @author John Bogovic <bogovicj@janelia.hhmi.org>
 * 
 */
public class LowerTriangularMatrixIterator implements Iterator< Long >
{

	protected int nrows; // width and height of matrix
	protected long numElements;

	protected long start; // starting point
	protected long n; // current index
	protected long count; // count
	protected int i, j; // current row and column in matrix

	public LowerTriangularMatrixIterator( int nrows, int start )
	{
		this.nrows = nrows;

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
		int[] tmp = idxToCoords( start );
		i = tmp[ 0 ];
		j = tmp[ 1 ];
	}

	@Override
	public boolean hasNext()
	{
		return (count < numElements);
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

	public int[] nextCoord()
	{
		Long n = next();
		return idxToCoords( n );
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
	public int coordsToIdx( int i, int j )
	{
		return i + nrows * j;
	}

	/**
	 * 
	 * @param n
	 *            the linear index
	 * @return array containing the { rowIndex, columnIndex }
	 */
	public int[] idxToCoords( long n )
	{
		int[] coords = new int[ 2 ];
		idxToCoords( n, coords );
		return coords;
	}

	/**
	 * 
	 * @param n
	 *            the linear index
	 * @return array containing the { rowIndex, columnIndex }
	 */
	public void idxToCoords( long n, int[] coords )
	{

		coords[ 0 ] = (int) (n % nrows);

		n -= coords[ 0 ];
		n /= nrows;

		coords[ 1 ] = (int) (n % nrows);
	}

	public static void main( String[] args )
	{
		LowerTriangularMatrixIterator ltmi = new LowerTriangularMatrixIterator(
				5 );

		System.out.println( "numel: " + ltmi.numElements );

		// System.out.println( "(0,1): " + utmi.coordsToIdx( 0, 1 ) );
		// System.out.println( "(1,0): " + utmi.coordsToIdx( 1, 0 ) );
		//
		// System.out.println( "(0): "
		// + ArrayUtil.printArray( utmi.idxToCoords( 0 ) ) );
		// System.out.println( "(1): "
		// + ArrayUtil.printArray( utmi.idxToCoords( 1 ) ) );
		// System.out.println( "(5): "
		// + ArrayUtil.printArray( utmi.idxToCoords( 5 ) ) );

		while ( ltmi.hasNext() )
		{
			System.out.println( ltmi.next() );
		}

		System.out.println( "done" );
	}
}
