package net.imglib2.algorithms.crack;

import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;

import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;

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
	
	private float eps = 0.01f;
	
	final protected int nx, ny;
	public int[] sz;
	public float getNx(){ return nx; }
	public float getNy(){ return ny; }
	
	final protected HashMap< AffineModel2D, ArrayList< PointMatch > > av = new HashMap< AffineModel2D, ArrayList< PointMatch > >();
	public HashMap< AffineModel2D, ArrayList< PointMatch > > getAV(){ return av; }
	final protected HashMap< PointMatch, ArrayList< AffineModel2D > > va = new HashMap< PointMatch, ArrayList< AffineModel2D > >();
	public HashMap< PointMatch, ArrayList< AffineModel2D > > getVA(){ return va; };

	final protected HashMap< AffineModel2D, Boolean > ac = new HashMap< AffineModel2D, Boolean >();
	public HashMap< AffineModel2D, Boolean > getAC(){ return ac; }
	
	final static protected PointFactory< Point > defaultPointFactory = new PointFactory< Point >()
	{
		final public Point createPoint( final float[] l )
		{
			return new Point( l );
		}
	};

	final static protected PointMatchFactory< PointMatch > defaultPointMatchFactory = new PointMatchFactory< PointMatch >()
	{
		
		final public PointMatch createPointMatch( final Point p1, final Point p2 )
		{
			return new PointMatch( p1, p2 );
		}

		
		final public PointMatch createPointMatch( final Point p1, final Point p2, final float w )
		{
			return new PointMatch( p1, p2, w );
		}
	};

	final static Logger logger = LogManager.getLogger( CrackTransformMesh.class.getName() );

	public CrackTransformMesh( int nx, int ny ){
		this.nx = nx;
		this.ny = ny;
		sz = new int[]{ nx, ny };
		
	}
	public CrackTransformMesh( int[] sz ){
		this.sz = sz;
		this.nx = sz[0];
		this.ny = sz[1];
	}
	
	public void reset(){
		av.clear();
		va.clear();
		ac.clear();
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
		
		float[][] ptsI   = null;
		float[][] ptsIp1 = null;
		
		for ( int i = 0; i < N - 1; i++ )
		{
			// save a bit of computation by reusing points
			// from the previous iteration
			if( i == 0){
				ptsI   = xfmPtList( ctrlPts[i], ctrlPtOffsets[i],     sz, eps ); 	// i 
				ptsIp1 = xfmPtList( ctrlPts[i+1], ctrlPtOffsets[i+1], sz, eps ); // i + 1	
			}else{
				ptsIp1 = xfmPtList( ctrlPts[i+1], ctrlPtOffsets[i+1], sz, eps ); // i + 1
			}
			
//			System.out.println(" ptsI:\n " + ArrayUtil.printArray( ptsI ));
//			System.out.println(" ");
//			System.out.println(" ptsIp1:\n " + ArrayUtil.printArray( ptsIp1 ));
			
			addTrianglesFromAdjacentPoints( ptsI, ptsIp1 );
			
			ptsI = ptsIp1; // (i+1) becomes (i) for the next iteration
		}
	}
	
	public static float[] intersection2d( float[] pt, float[] off, int[] sz, boolean neg )
	{
		float[] Xpt = new float[pt.length];
		
		// projection to 4 faces
		float[][] faceVecs = createVectorsToFaces( pt, sz );
		float maxDotR = -1f;
		int k = -1;
		for ( int i=0; i<4; i++)
		{
			float dot = off[0] * faceVecs[i][0] + off[1] * faceVecs[i][1];
			
			if( neg ){ dot *= -1; }
			
			if( dot <= 0 ){
				continue;
			}
			double mag = Math.sqrt( (double)ArrayUtil.sumSquares( faceVecs[i] ) );
			double dotRatio = dot/mag;
			
			if( dotRatio > maxDotR ){
				maxDotR = (float)dotRatio;
				k = i;
			}
		}
		
		float r = 0; 
		if ( faceVecs[k][0] == 0f )
		{
			r = (faceVecs[k][1] / maxDotR);
		}else{
			r = (faceVecs[k][0] / maxDotR);
		}
		
		Xpt[0] = pt[0] + r * off[0];
		Xpt[1] = pt[1] + r * off[1];
		
		return Xpt;
	}
	
	public static float[][] createVectorsToFaces( float[] pt, int[] sz ){
		float[][] vecs = new float[4][2];
		
		//v1
		vecs[0][0] = - pt[0];
		//v2
		vecs[1][1] = - pt[1];
		//v3
		vecs[2][0] = sz[0] - pt[0];
		//v4
		vecs[3][1] = sz[1] - pt[1];
		return vecs;
	}

	public static float[][] xfmPtList( float[] pt, float[] offset, int[] sz, float eps )
	{
		float[][] ptList = new float[12][];
		
		ptList[0] = intersection2d( pt, offset, sz, true );
		ptList[1] = linCombo( 1f, ptList[0], -1f-eps, offset );
		
		ptList[2] = linCombo( 1f, pt,  -2*eps, offset );
		ptList[3] = linCombo( 1f, pt, -1f-eps, offset );
		
		ptList[4] = linCombo( 1f, pt, -eps, offset );
		ptList[5] = linCombo( 1f, pt,  -1f, offset );
		
		ptList[6] = linCombo( 1f, pt, eps, offset );
		ptList[7] = linCombo( 1f, pt,  1f, offset );
		
		ptList[8] = linCombo( 1f, pt,  2*eps, offset );
		ptList[9] = linCombo( 1f, pt, 1f+eps, offset );
		
		ptList[10] = intersection2d( pt, offset, sz, false );
		ptList[11] = linCombo( 1f, ptList[10], 1f+eps, offset );
		
		return ptList;
	}
	
	public static float[][] xfmPtListDeform( float[] pt, float[] offset, int[] sz, float eps )
	{
		float[][] ptList = new float[12][];
		
		ptList[0] = intersection2d( pt, offset, sz, true );
		ptList[1] = ptList[0];
		
		ptList[2] = linCombo( 1f, pt,  -2*eps, offset );
		ptList[3] = linCombo( 1f, pt, -1f-eps, offset );
		
		ptList[4] = linCombo( 1f, pt, -eps, offset );
		ptList[5] = linCombo( 1f, pt,  -1f, offset );
		
		ptList[6] = linCombo( 1f, pt, eps, offset );
		ptList[7] = linCombo( 1f, pt,  1f, offset );
		
		ptList[8] = linCombo( 1f, pt,  2*eps, offset );
		ptList[9] = linCombo( 1f, pt, 1f+eps, offset );
		
		ptList[10] = intersection2d( pt, offset, sz, false );
		ptList[11] = ptList[10];
		
		return ptList;
	}
	
	public void addTrianglesFromAdjacentPoints( float[][] ptList1, float[][] ptList2 )
	{
		if( ptList1.length != ptList2.length ){
			logger.error("Point lists must be the same length.");
			return;
		}

		int N = ptList1.length;
		ArrayList<PointMatch> triangle = null;

		// local points
		float[] p1 = null;
		float[] p2 = null;
		float[] p3 = null;

		// matching world points
		float[] m1 = null;
		float[] m2 = null;
		float[] m3 = null;

		// there are N-1 spaces between N points
		// and two triangles per space
		boolean even = true;
		for( int i=0; i<N-2; i+=2 )
		{
			// Illustration showing first and second triangles in the loop
			// T1 is the 'odd'  triangle
			// T2 is the 'even' triangle,  
			//
			//  ptList[i+2] #-------# ptList2[i+2]   ptList[i+3] #-------# ptList2[i+3] 
			//              |      /|                            |      /|
			//              | T1  / |                            | T1  / |
			//   [local]    |    /  |         -------->          |    /  |    [world]
			//              |   /   |                            |   /   |
			//              |  /    |                            |  /    |
			//              | / T2  |                            | / T2  |
			//    ptList[i] #-------# ptList2[i]     ptList[i+2] #-------# ptList2[i+2] 

			int j = i + 1;

			// two triangles per i
			for (int k = 0; k<2; k++)
			{	
				triangle = new ArrayList<PointMatch>(3);
				
				if( even ){  
					logger.debug("even triangle");
					// even triangles use two points from ptList2
					p1 = ptList1[i];	
					p2 = ptList2[i];	
					p3 = ptList2[i+2];	

					m1 = ptList1[j];	
					m2 = ptList2[j];	
					m3 = ptList2[j+2];	

				}else{
					logger.debug("odd  triangle");
					// odd triangles use two points from ptList1
					p1 = ptList1[i];	
					p2 = ptList1[i+2];	
					p3 = ptList2[i+2];	 

					m1 = ptList1[j];	
					m2 = ptList1[j+2];	
					m3 = ptList2[j+2];	 
				}

				triangle.add( 
						new PointMatch(
								new Point( p1 ),
								new Point( p1, m1 )
						));

				triangle.add( 
						new PointMatch(
								new Point( p2 ),
								new Point( p2, m2 )
						));

				triangle.add( 
						new PointMatch(
								new Point( p3 ),
								new Point( p3, m3 )
						));

				if ( i == (N/2) - 2 ){
					addTriangle( triangle, true );
				}else{
					addTriangle( triangle );
				}
				
				
//				System.out.println("\n added triangle:\nLocal:\n  " +
//						triangle.get(0).getP1().getL()[0] + " " + triangle.get(0).getP1().getL()[1] + "\n  " +
//						triangle.get(1).getP1().getL()[0] + " " + triangle.get(1).getP1().getL()[1] + "\n  " +
//						triangle.get(2).getP1().getL()[0] + " " + triangle.get(2).getP1().getL()[1] + "\n" +
//						"World:\n  " +
//						triangle.get(0).getP2().getW()[0] + " " + triangle.get(0).getP2().getW()[1] + "\n  " +
//						triangle.get(1).getP2().getW()[0] + " " + triangle.get(1).getP2().getW()[1] + "\n  " +
//						triangle.get(2).getP2().getW()[0] + " " + triangle.get(2).getP2().getW()[1] + "\n"
//					);

				// switch from even to odd + vice versa
				even = !even;
			}


		}
	}
	

	public PointMatch makeMatch( float[] pt, float[] off, float m1, float c1, float m2, float c2 ){
		
		float[] ptlc = ArrayUtil.clone(pt);
		linComboInPlace( 1, ptlc, m1, off, c1 );
		
		float[] mtch = ArrayUtil.clone(pt);
		linComboInPlace( 1, mtch, m2, off, c2 );
		
		return new PointMatch(
				new Point( ptlc ),
				new Point( ptlc, mtch )
				);
	}
	
	public PointMatch makeMatch( float[] pt, float[] off, float amt0, float amt1, float eps ){
		
		float[] ptlc = ArrayUtil.clone(pt);
		linComboInPlace( 1, ptlc, amt0, off );
		
		float[] mtch = ArrayUtil.clone(pt);
		linComboInPlace( 1, mtch, amt1, off, eps );
		
		return new PointMatch(
				new Point( ptlc ),
				new Point( ptlc, mtch )
				);
	}

	public PointMatch makeMatch( float[] pt, float[] off, float amt0, float amt1 ){
		
		float[] ptlc = ArrayUtil.clone(pt);
		linComboInPlace( 1, ptlc, amt0, off );
		
		float[] mtch = ArrayUtil.clone(pt);
		linComboInPlace( 1, mtch, amt1, off );
		
		return new PointMatch(
				new Point( ptlc ),
				new Point( ptlc, mtch )
				);
	}

	/**
	 * returns out = a*v1 = b*v2
	 * @param a
	 * @param v1
	 * @param b
	 * @param v2
	 */
	private static float[] linCombo( float a, float[] v1, float b, float[] v2 ){
		int N = v1.length;
		float[] out = v1.clone();
		linComboInPlace(a, out, b, v2);
		return out;
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
	
	private static float[] linCombo( float a, float[] v1, float b, float[] v2, float c ){
		float[] out = v1.clone();
		linComboInPlace( a, out, b, v2, c);
		return out;
	}
	
	/**
	 * sets v1 = a*v1 = b*v2 + c
	 * @param a
	 * @param v1
	 * @param b
	 * @param v2
	 */
	private static void linComboInPlace( float a, float[] v1, float b, float[] v2, float c ){
		int N = v1.length;
		for( int i = 0; i < N; i++ ){
			v1[i] = a * v1[i] + b * v2[i] + c;
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
	public void addTriangle( final ArrayList< PointMatch > t , boolean inCrack )
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
			if ( !va.containsKey( pm ) ){
				va.put( pm, new ArrayList< AffineModel2D >() );
			}
			va.get( pm ).add( m );
		}
		
		ac.put( m, inCrack );
	}	
	
	public void addTriangle( final ArrayList< PointMatch > t ){
		addTriangle( t, false );
	}

	public float[] apply(float[] location) {
		assert location.length == 2 : "2d transform meshs can be applied to 2d points only.";
		
		final float[] transformed = location.clone();
		applyInPlace( transformed );
		return transformed;
	}

	public void applyInPlace(float[] location) {
		assert location.length == 2 : "2d transform meshs can be applied to 2d points only.";
		
		final Set< AffineModel2D > s = av.keySet();
		for ( final AffineModel2D ai : s )
		{
			final ArrayList< PointMatch > pm = av.get( ai );
			if ( TransformMesh.isInSourcePolygon( pm, location ) )
			{
				logger.info(" in polygon: ");
				ai.applyInPlace( location );
				return;
			}
		}
		// else return the identity 
	}

	public float[] applyInverse(float[] point)
			throws NoninvertibleModelException {
		assert point.length == 2 : "2d transform meshs can be applied to 2d points only.";
		
		final float[] transformed = point.clone();
		applyInverseInPlace( transformed );
		return transformed;
	}

	public void applyInverseInPlace(float[] point)
			throws NoninvertibleModelException {
		assert point.length == 2 : "2d transform meshs can be applied to 2d points only.";
		
		final Set< AffineModel2D > s = av.keySet();
		for ( final AffineModel2D ai : s )
		{
			final ArrayList< PointMatch > pm = av.get( ai );
			if ( TransformMesh.isInConvexTargetPolygon( pm, point ) )
			{
				ai.applyInverseInPlace( point );
				return;
			}
		}
//		throw new NoninvertibleModelException( "Noninvertible location ( " + point[ 0 ] + ", " + point[ 1 ] + " )" );
	}

	public InvertibleCoordinateTransform createInverse() {
		// TODO Auto-generated method stub
		return null;
	}
	
	public String toString(){
		String out = "CrackTransformMesh with:\n";
		Set<Entry<AffineModel2D, ArrayList<PointMatch>>> triangles = av.entrySet();
		
		
		out += triangles.size() + " triangles\n";
		Iterator<Entry<AffineModel2D, ArrayList<PointMatch>>> triangleIt = triangles.iterator();
		
		while( triangleIt.hasNext() )
		{
			Entry<AffineModel2D, ArrayList<PointMatch>> ent = triangleIt.next();
			Boolean inCrack = ac.get( ent.getKey() );
			if( inCrack.booleanValue() ){
				out +=  " CRACK \n";
			}else{
				out +=  " NO CRACK \n";
			}
			out += 	    "  " + ArrayUtil.printArray(ent.getValue().get(0).getP1().getL())
					+ " -> " + ArrayUtil.printArray(ent.getValue().get(0).getP2().getW()) + "\n";
			out += 	    "  " + ArrayUtil.printArray(ent.getValue().get(1).getP1().getL())
					+ " -> " + ArrayUtil.printArray(ent.getValue().get(1).getP2().getW()) + "\n";
			out += 	    "  " + ArrayUtil.printArray(ent.getValue().get(2).getP1().getL())
					+ " -> " + ArrayUtil.printArray(ent.getValue().get(2).getP2().getW()) + "\n";
			out += " *** \n";
		}
		return out;
	}
	public String toStringMatlab(){
		String out = "CrackTransformMesh with:\n";
		Set<Entry<AffineModel2D, ArrayList<PointMatch>>> triangles = av.entrySet();
		out += triangles.size() + " triangles\n";
		Iterator<Entry<AffineModel2D, ArrayList<PointMatch>>> triangleIt = triangles.iterator();
		
		while( triangleIt.hasNext() )
		{
			Entry<AffineModel2D, ArrayList<PointMatch>> ent = triangleIt.next();
			Boolean inCrack = ac.get( ent.getKey() );
			if( inCrack.booleanValue() ){
				out +=  " CRACK \n";
			}else{
				out +=  " NO CRACK \n";
			}
			out += 	    "  " + ArrayUtil.printArray(ent.getValue().get(0).getP1().getL())
					+ " " + ArrayUtil.printArray(ent.getValue().get(0).getP2().getW()) + ";...\n";
			out += 	    "  " + ArrayUtil.printArray(ent.getValue().get(1).getP1().getL())
					+ " " + ArrayUtil.printArray(ent.getValue().get(1).getP2().getW()) + ";...\n";
			out += 	    "  " + ArrayUtil.printArray(ent.getValue().get(2).getP1().getL())
					+ " " + ArrayUtil.printArray(ent.getValue().get(2).getP2().getW()) + ";...\n";
		}
		return out;
	}
	
	public static void linComboTest(){
		float[] v1 = new float[]{ 1,2,3};
		float[] v2 = new float[]{ 1,1,1};
		linComboInPlace( 1f, v1, -0.5f, v2);
		System.out.println( " " + ArrayUtil.printArray( v1 ));
	}
	
	
	public static void resampTest(){
//		String fn = "/Users/bogovicj/Documents/learning/advanced-imglib2/images/bee-1.tif";
//		Img<FloatType> im = ImagePlusAdapter.convertFloat( IJ.openImage(fn) );
		
		String fnOut = "/Users/bogovicj/Documents/projects/crackSim/grad1/grad_meshCrack_push.tif";
		
		int[] sz = new int[]{ 200, 200 };
		Img<FloatType> im = ImgOps.createGradientImgY( sz, new FloatType());
		
		ImageProcessor ipin = new FloatProcessor( sz[0], sz[1]);
		ImageProcessor ipout = new FloatProcessor( sz[0], sz[1]);
		try {
			ImgOps.copyToImageProcessor2dFloat(im, ipin);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
//		float[][] pts = new float[][]{ {100,100},  {105,100}, {145,100}, {150,100}};
//		float[][] off = new float[][]{ {0f,0.01f}, {0f,10f},  {0f,11f},  {0f,0.01f}};
		
//		float[][] pts = new float[][]{ {100,100},  {125,100}, {150,100}};
//		float[][] off = new float[][]{ {0f,0.01f}, {0f,10f},  {0f,0.01f}};
		
		float[][] pts = new float[][]{ {100,100}, {115,100}, {135,100}, {150,100} };
		float[][] off = new float[][]{ {0f,0.2f}, {0f,6f} ,   {0f,5f}, {0f,0.2f}  };
		
		CrackTransformMesh crackMesh = new CrackTransformMesh( 200, 200 );
		crackMesh.fromCrackParam( pts , off );
		
//		System.out.println( " crack mesh :\n" + crackMesh.toString() );
		
//		Set<Entry<AffineModel2D, Boolean>> ents = crackMesh.getAC().entrySet();
//		Iterator<Entry<AffineModel2D, Boolean>> entIt = ents.iterator();
//		while ( entIt.hasNext() ){
//			Entry<AffineModel2D, Boolean> e = entIt.next();
//			if(e.getValue().booleanValue()){
//				System.out.println(" an xfm in a crack! ");
//			}else{
//				System.out.println(" non-crack xfm ");
//			}
//		}
		
		CrackTransformMeshMapping<CrackTransformMesh> map 
			= new CrackTransformMeshMapping<CrackTransformMesh>( crackMesh );
		
		map.mapInterpolated( ipin, ipout, 1 );		
		ImagePlus impout = new ImagePlus( "hi", ipout );
		IJ.save( impout, fnOut);
		
	}
	
	public static void mapTest() throws NoninvertibleModelException{
		float[][] pts = new float[][]{ {100,100}, {150,100}};
		float[][] off = new float[][]{ {0f,10f}, {0f,10f}};
		
		CrackTransformMesh crackMesh = new CrackTransformMesh( 256, 256 );
		crackMesh.fromCrackParam( pts , off );
		
		System.out.println( "" + crackMesh );
		
		CrackTransformMeshMapping<CrackTransformMesh> map 
			= new CrackTransformMeshMapping<CrackTransformMesh>( crackMesh );
	
		float[] pt = new float[]{125f, 105f}; 
//		System.out.println(" mapped " + ArrayUtil.printArray(pt) + " to " +
//				 ArrayUtil.printArray(crackMesh.apply(pt)) );
//		
//		pt = new float[]{50f, 50f}; 
//		System.out.println(" mapped " + ArrayUtil.printArray(pt) + " to " +
//				 ArrayUtil.printArray(crackMesh.apply(pt)) );
		
		
		System.out.println( "\n*\n*\n" );
		for ( float y = 85f; y<115f; y++ ){
			pt = new float[]{125f, y};
			float[] ptm = crackMesh.apply(pt);
			System.out.println(" map" + pt[1] + " to " + ptm[1] );
		}
		
	}
	public static void debugTriangulation(){
		float[][] pts = new float[][]{ {100,100}, {150,100}};
		float[][] off = new float[][]{ {0f,10f},  {0f,10f}};
		
		CrackTransformMesh crackMesh = new CrackTransformMesh( 200, 200 );
		crackMesh.fromCrackParam( pts , off );
		
		System.out.println( "" + crackMesh );
		
		Set<Entry<AffineModel2D, ArrayList<PointMatch>>> entries = crackMesh.getAV().entrySet();
		Iterator<Entry<AffineModel2D, ArrayList<PointMatch>>> it = entries.iterator();
		System.out.println(" \n"  + "Local:");
		while ( it.hasNext() ){
			ArrayList<PointMatch> pm = it.next().getValue();
			
			float[] a = pm.get( 0 ).getP1().getL();
			float ax = a[ 0 ];
			float ay = a[ 1 ];
			float[] b = pm.get( 1 ).getP1().getL();
			float bx = b[ 0 ];
			float by = b[ 1 ];
			float[] c = pm.get( 2 ).getP1().getL();
			float cx = c[ 0 ];
			float cy = c[ 1 ];
			
			System.out.println(" "  +
					""+ax + "," + ay + "," + bx + ","+by+"," + cx+","+cy +";...");
		}
		System.out.println(" \n"  + "World:");
		it = entries.iterator();
		while ( it.hasNext() ){
			ArrayList<PointMatch> pm = it.next().getValue();
			
			float[] a = pm.get( 0 ).getP2().getW();
			float ax = a[ 0 ];
			float ay = a[ 1 ];
			float[] b = pm.get( 1 ).getP2().getW();
			float bx = b[ 0 ];
			float by = b[ 1 ];
			float[] c = pm.get( 2 ).getP2().getW();
			float cx = c[ 0 ];
			float cy = c[ 1 ];
			
			System.out.println(" "  +
					""+ax + "," + ay + "," + bx + ","+by+"," + cx+","+cy +";...");
		}
	}
	

	
	public static void testFaceDistance(){
		
		int[]   sz = new int[]  { 15, 15 };
		float[] pt = new float[]{ 12, 12 };
		float[] offset = ArrayUtil.normalizeLength( new float[]{ 0.4f, 0.7f} );
		
//		float[][] faceVecs = createVectorsToFaces( pt, sz);
//		System.out.println(" " + ArrayUtil.printArray(faceVecs));
//		System.out.println(" \n*\n*\n " );
		
		float[] Xpt = intersection2d( pt, offset, sz, false );
		System.out.println("\n pt:     " + ArrayUtil.printArray( pt ));
		System.out.println(" offset: " + ArrayUtil.printArray( offset ));
		System.out.println(" Xpt " + ArrayUtil.printArray(Xpt));
		
		
		offset = ArrayUtil.normalizeLength( new float[]{ 0.7f, 0.4f} );
		Xpt = intersection2d( pt, offset, sz, false );
		System.out.println("\n pt:     " + ArrayUtil.printArray( pt ));
		System.out.println(" offset: " + ArrayUtil.printArray( offset ));
		System.out.println(" Xpt " + ArrayUtil.printArray(Xpt));
		System.out.println(""); 
		
	}
	
	public static void ptListTest(){
		float[] pt = new float[]{100,100};	
		float[] pt2 = new float[]{110,100};
		
		float[] offset = new float[]{0f, 10f};
		int[] sz = new int[]{ 256, 256};
		float eps = 0.01f;
		
		float[][] ptList  = xfmPtList( pt, offset, sz, eps );
		float[][] ptList2 = xfmPtList( pt2, offset, sz, eps );
		
		System.out.println( ArrayUtil.printArray( ptList ));
		System.out.println( "" );
		System.out.println( ArrayUtil.printArray( ptList2 ));

//		float[][] pts = new float[][]{ {100,100}, {150,100}};
//		float[][] off = new float[][]{ {0f,10f}, {0f,10f}};
//		
//		CrackTransformMesh crackMesh = new CrackTransformMesh( 200, 200, 10f, 10f );
//		crackMesh.fromCrackParam( pts , off );
//		
//		System.out.println("\n" + crackMesh.toStringMatlab() );
		
	}
	
	public static void main(String[] args){

//		ptListTest();
		
//		debugTriangulation();
		
//		testFaceDistance();
		
		resampTest();
		
//		try {
//			mapTest();
//		} catch (NoninvertibleModelException e) {
//			e.printStackTrace();
//		}
	
		System.out.println("Finished");
		System.exit(0);
	}
	

}
