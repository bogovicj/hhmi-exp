package net.imglib2.algorithms.crack;

import java.util.ArrayList;
import java.util.HashMap;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import edu.jhu.ece.iacl.utility.ArrayUtil;
import mpicbg.models.*;

/**
 * Inspired by {@link TransformMesh}, but specific to crack.
 * 
 * 
 * 
 * @author John Bogovic
 * @author Stephan Saalfeld
 */
public class CrackTransformMesh implements InvertibleCoordinateTransform 
{
	
	private float eps = 0.3f;
	
	final protected int nx, ny;
	public float getNx(){ return nx; }
	public float getNy(){ return ny; }
	
	final protected float width, height;
	public float getWidth(){ return width; }
	public float getHeight(){ return height; }
	
	final protected HashMap< AffineModel2D, ArrayList< PointMatch > > av = new HashMap< AffineModel2D, ArrayList< PointMatch > >();
	public HashMap< AffineModel2D, ArrayList< PointMatch > > getAV(){ return av; }
	final protected HashMap< PointMatch, ArrayList< AffineModel2D > > va = new HashMap< PointMatch, ArrayList< AffineModel2D > >();
	public HashMap< PointMatch, ArrayList< AffineModel2D > > getVA(){ return va; };

	final static protected PointFactory< Point > defaultPointFactory = new PointFactory< Point >()
	{
		@Override
		final public Point createPoint( final float[] l )
		{
			return new Point( l );
		}
	};

	final static protected PointMatchFactory< PointMatch > defaultPointMatchFactory = new PointMatchFactory< PointMatch >()
	{
		@Override
		final public PointMatch createPointMatch( final Point p1, final Point p2 )
		{
			return new PointMatch( p1, p2 );
		}

		@Override
		final public PointMatch createPointMatch( final Point p1, final Point p2, final float w )
		{
			return new PointMatch( p1, p2, w );
		}
	};

	final static Logger logger = LogManager.getLogger( CrackTransformMesh.class.getName() );

	public CrackTransformMesh( int nx, int ny, float width, float height){
		this.nx = nx;
		this.ny = ny;
		
		this.width  = width;
		this.height = height;
	}
	
	public void reset(){
		av.clear();
		va.clear();
	}
	
	public void fromCrackParam( float[][] ctrlPts, float[][] ctrlPtOffsets )
	{
		int N = ctrlPts.length;
		// validate inputs
		if ( ctrlPtOffsets.length != N ){
			logger.error( "Inputs must be the same length" );
			return;
		}
		if ( ctrlPts[0].length != 2 || ctrlPtOffsets[0].length != 2 ){
			logger.error( "Inputs must be 2d points" );
			return;
		}
		
		float[] curPt = null;
		float[] nxtPt = null;
		
		for ( int i = 0; i < N - 1; i++ )
		{
			
		}
		
	}
	
	public void add( float[] curPt, float[] nxtPt ){
		float[] p1 = new float[2];
		float[] p2 = new float[2];
		float[] p3 = new float[2];
		
		
		
	}
	
	public PointMatch makeMatch( float[] pt, float[] off, float amt ){
		
		float[] mtch = ArrayUtil.clone(pt);
		linComboInPlace( 1, pt, amt, off );
		
		return new PointMatch(
				new Point( pt ),
				new Point( mtch )
				);
	}

	/**
	 * sets v1 = a*v1 = b*v2
	 * @param a
	 * @param v1
	 * @param b
	 * @param v2
	 */
	private static void linComboInPlace( float a, float[] v1, float b, float[] v2 ){
		int N = v1.length;
		for( int i = 0; i < N; i++ ){
			v1[i] = a * v1[i] + b * v2[i];
		}
	}
	
	/**
	 * Add a triangle defined by 3 PointMatches that defines an
	 * AffineTransform2D.
	 * 
	 * @param t
	 *            3 PointMatches (will not be copied, so do not reuse this
	 *            list!)
	 */
	public void addTriangle( final ArrayList< PointMatch > t )
	{
		final AffineModel2D m = new AffineModel2D();
		try
		{
			m.fit( t );
		}
		catch ( final NotEnoughDataPointsException e ) { e.printStackTrace(); }
		catch ( final IllDefinedDataPointsException e ) { e.printStackTrace(); }
		av.put( m, t );
		
		for ( final PointMatch pm : t )
		{
			if ( !va.containsKey( pm ) )
				va.put( pm, new ArrayList< AffineModel2D >() );
			va.get( pm ).add( m );
		}
	}	
	
	@Override
	public float[] apply(float[] location) {
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public void applyInPlace(float[] location) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public float[] applyInverse(float[] point)
			throws NoninvertibleModelException {
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public void applyInverseInPlace(float[] point)
			throws NoninvertibleModelException {
		// TODO Auto-generated method stub
		
	}
	@Override
	public InvertibleCoordinateTransform createInverse() {
		// TODO Auto-generated method stub
		return null;
	}
	
	
	public static void main(String[] args){
		float[] v1 = new float[]{ 1,2,3};
		float[] v2 = new float[]{ 1,1,1};
		linComboInPlace( 1f, v1, -0.5f, v2);
		System.out.println( " " + ArrayUtil.printArray( v1 ));
	}
	

}
