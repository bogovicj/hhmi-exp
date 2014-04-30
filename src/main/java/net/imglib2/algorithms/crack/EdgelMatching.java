package net.imglib2.algorithms.crack;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithms.edge.EdgelTools;
import net.imglib2.collection.KDTree;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.imageplus.ImagePlusImgFactory;
import net.imglib2.neighborsearch.KNearestNeighborSearchOnKDTree;
import net.imglib2.neighborsearch.RadiusNeighborSearchOnKDTree;
import net.imglib2.realtransform.InverseRealTransform;
import net.imglib2.realtransform.RealTransformRandomAccessible;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;
import net.imglib2.util.LinAlgHelpers;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import edu.jhu.ece.iacl.utility.ArrayUtil;

/**
 *
 */
public class EdgelMatching<T extends NativeType<T> & RealType<T>>{

	Logger logger = LogManager.getLogger(EdgelMatching.class.getName());

	public enum SearchTypes { RADIUS, COUNT };
	
	private Img<T> img;
	private Img<T> mask;
	private int[] patchSize;

	private int numEdgels;
	private ArrayList<Edgel> 	edgels;
	private KDTree<Edgel>		edgelTree;

	private HashMap<EdgelPair, Double> edgelAffinities;

	// search parameters 
	private double edgelSearchRadius = 5;
	private int    edgelSearchCount  = 15;
	private SearchTypes	search = SearchTypes.RADIUS;
	private KNearestNeighborSearchOnKDTree<Edgel>  countSearch;
	private RadiusNeighborSearchOnKDTree<Edgel>    radiusSearch;
	
	int i, j;
	
	// comparison parameters
	private double lam = 0.5;
	
	// temporary storage
	Img<T> depth1;
	Img<T> depth2;

	public EdgelMatching( Img<T> img, Img<T> mask, int[] patchSize)
	{
		this.img = img;
		this.mask = mask;

		this.patchSize = patchSize;
	}
	
	public void setSearchType( SearchTypes searchType ){
		this.search = searchType;
	}

	public double getEdgelSearchRadius() {
		return edgelSearchRadius;
	}

	public void setEdgelSearchRadius(double edgelSearchRadius) {
		this.edgelSearchRadius = edgelSearchRadius;
	}

	public int getEdgelSearchCount() {
		return edgelSearchCount;
	}

	public void setEdgelSearchCount(int edgelSearchCount) {
		this.edgelSearchCount = edgelSearchCount;
	}

	/**
	 * @return the edgels
	 */
	public ArrayList<Edgel> getEdgels() {
		return edgels;
	}

	/**
	 * @param edgels the edgels
	 */
	public void setEdgels(ArrayList<Edgel> edgels) 
	{
		this.edgels = edgels;
		this.numEdgels = edgels.size();
		
	}

	/**
	 * Returns edgels that are likely to match the i^th edgel based 
	 * on the local crack geometry.
	 */
	public ArrayList<Edgel> candidateEdgels(Edgel e)
	{
		switch ( search )
		{
		case COUNT:
			return candidateEdgelsK(e);
			
		default:
			return candidateEdgelsRadius(e);
		}
		
	}

	public ArrayList<Edgel> candidateEdgelsRadius(Edgel e)
	{
		ArrayList<Edgel> out = new ArrayList<Edgel>();
		
		radiusSearch.search( e, getEdgelSearchRadius(), false);
		for( int i=0; i<radiusSearch.numNeighbors(); i++ )
		{
			out.add( radiusSearch.getSampler(i).get() );
		}
		
		return out;
	}
	
	public ArrayList<Edgel> candidateEdgelsK(Edgel e)
	{
		ArrayList<Edgel> out = new ArrayList<Edgel>( getEdgelSearchCount() );
		
		countSearch.search( e );
		for( int i=0; i<countSearch.getK(); i++ )
		{
			out.add( countSearch.getSampler(i).get() );
		}
		
		return out;
	}
	
	/**
	 * Removes an edgel from the candidate list if its normal direction 
	 * points in the opposite direction as the reference Edgel e. 
	 * ( i.e. if scalar product with ref is negative ).  
	 * @param e
	 * @param candidates
	 */
	public void filterEdgelsByNormal(Edgel e, List<Edgel> candidates)
	{
		
		for(int i=candidates.size()-1; i>=0; i--)
		{
			double dot = LinAlgHelpers.dot(e.getGradient(), candidates.get(i).getGradient());
			logger.trace(" dot " + dot);
			if(  dot >= 0 )
				candidates.remove(i);
		}
	}


	/**
	 * Computes a likelihood that edgels i and j are a match.
	 */
	public double edgelAffinities( Edgel e1, Edgel e2 )
	{
		
		RealTransformRandomAccessible<T, InverseRealTransform> e1View = 
				EdgelTools.edgelToView(e1, img, patchSize);
		
		RealTransformRandomAccessible<T, InverseRealTransform> e2View = 
				EdgelTools.edgelToView(e2, img, patchSize);
		
		CrackCorrection.computeCrackDepthNormalMask(e1, mask, patchSize, depth1);
		CrackCorrection.computeCrackDepthNormalMask(e2, mask, patchSize, depth2);
		
		double depthSSD = 0;
		
//		Cursor<FloatType> curs = depth1.cursor();
//		while( curs.hasNext() )
//		{
//			curs.fwd();
//			
//			
//		}
		
		return -1;
	}
	
	/**
	 * Computes an affinity between two edgels using only geometrical information.
	 * <P>
	 * dist(e1,e1) + lam * dot
	 * @param e1 first edgel
	 * @param e2 second edgel
	 * @param lam contribution of normal vector
	 * @return
	 */
	
	public double edgelAffinitiesGeom ( Edgel e1, Edgel e2, double lam )
	{
		double dot = -LinAlgHelpers.dot( e1.getGradient(), e2.getGradient() );
		
		double[] p1 = new double[e1.numDimensions()];
		e1.localize(p1);
		double[] p2 = new double[e1.numDimensions()];
		e2.localize(p2);
		
		double aff = LinAlgHelpers.squareLength( 
				ArrayUtil.subtract( p1, p2 ) );
		
		return ( 1 - lam ) * aff + ( lam * aff ) * dot ;
	}

	/**
	 *
	 */
	public void testAffinitiesReg()
	{
		logger.info("Computing edgel affinities");
		
		ArrayImgFactory<T> factory = new ArrayImgFactory<T>(); 
		depth1 = factory.create(ArrayUtil.toLong(patchSize), img.firstElement());
		depth2 = factory.create(ArrayUtil.toLong(patchSize), img.firstElement());
		
		edgelTree = new KDTree<Edgel>( edgels, edgels );
		switch ( search )
		{
		case COUNT:
			countSearch = new  KNearestNeighborSearchOnKDTree<Edgel>( edgelTree, getEdgelSearchCount());
			logger.info("" + getEdgelSearchCount() +"-Nearest neighbor search");
		default:
			radiusSearch = new  RadiusNeighborSearchOnKDTree<Edgel>( edgelTree );
			logger.info("" + getEdgelSearchRadius()+  " radius search");
		}
		
		edgelAffinities = new HashMap<EdgelPair,Double>();

//		int i = 12725;
		
		
//		i = cc.edgelIdxNearest(new double[]{67,290,13});
		
		Edgel e = edgels.get(i);
		logger.debug("i: " + i + "   " + e);

		ArrayList<Edgel> candidateEdgels = candidateEdgels(e);
		logger.debug(" " + candidateEdgels.size() + " matches after search");

		filterEdgelsByNormal(e, candidateEdgels);

		logger.debug(" " + candidateEdgels.size() + " matches after filtering.");

		if( candidateEdgels.size() == 0 ){ 
			logger.error(" no edgel matches left ");
			return;
		}
		
	}
	

	/**
	 *
	 */
	public void computeAllAffinities()
	{
		logger.info("Computing edgel affinities");
		
		ArrayImgFactory<T> factory = new ArrayImgFactory<T>(); 
		depth1 = factory.create(ArrayUtil.toLong(patchSize), img.firstElement());
		depth2 = factory.create(ArrayUtil.toLong(patchSize), img.firstElement());
		
		edgelTree = new KDTree<Edgel>( edgels, edgels );
		switch ( search )
		{
		case COUNT:
			countSearch = new  KNearestNeighborSearchOnKDTree<Edgel>( edgelTree, getEdgelSearchCount());
			logger.info("" + getEdgelSearchCount() +"-Nearest neighbor search");
		default:
			radiusSearch = new  RadiusNeighborSearchOnKDTree<Edgel>( edgelTree );
			logger.info("" + getEdgelSearchRadius()+  " radius search");
		}
		
		edgelAffinities = new HashMap<EdgelPair,Double>();

		int i = 0;
		for (Edgel e : edgels )
		{
			logger.debug("edgel " + i + " of " + numEdgels);
			
			ArrayList<Edgel> candidateEdgels = candidateEdgels(e);
			logger.debug(" " + candidateEdgels.size() + " matches after search");
			
			filterEdgelsByNormal(e, candidateEdgels);
			
			logger.debug(" " + candidateEdgels.size() + " matches after filtering.");

//			tabulateAffinities( e, candidateEdgels );
			
			int j = maxAffinityEdgelIdx( e, candidateEdgels );
			debugVisEdgelView( i, e, candidateEdgels.get(j) );
			
			
//			debugVisEdgelMatches( i, e, candidateEdgels );
			
			i++;
			
			// for debug 
			if( i > 20 ) break;
			
		}

	}
	
	public int maxAffinityEdgelIdx( Edgel e, List<Edgel> matches)
	{
		double maxAffinity = Double.MIN_VALUE;
		int    idxOut = -1;
		
		for( int i = 0 ; i<matches.size(); i++)
		{
			Edgel ec = matches.get(i);
			double affinity = edgelAffinitiesGeom( e, ec, 0.5 ); 
			if ( affinity > maxAffinity ){
				maxAffinity = affinity;
				idxOut = i;
			}
		}
		return idxOut;
	}
	
	public void tabulateAffinities( Edgel e, List<Edgel> matches)
	{
		for( Edgel ec : matches )
		{
			EdgelPair pair = new EdgelPair(e, ec);
			
			// if we already have this affinity skip the computation
			if( edgelAffinities.containsKey(pair) ) continue;	

			double affinity = edgelAffinities( e, ec ); 

			edgelAffinities.put(pair, affinity);
		}
	}
	
	
	public void debugVisEdgelView(int i, Edgel e, Edgel f){
		RandomAccessible<T> ev = 
				EdgelTools.edgelToView(e, img, patchSize);
		
		RandomAccessible<T> fv = 
				EdgelTools.edgelToView(f, img, patchSize);
		
//		ArrayImgFactory<T> ubfactory = new ArrayImgFactory<T>();
		ImagePlusImgFactory<T> ipfactory = new ImagePlusImgFactory<T>(); 
		
		Img<T> ePatchImg = ipfactory.create(patchSize, mask.firstElement());
		Img<T> fPatchImg = ipfactory.create(patchSize, mask.firstElement());
		
		logger.debug( " copying " );
		ImgOps.copyInto(ev, ePatchImg);
		ImgOps.copyInto(fv, fPatchImg);
		logger.debug( " done copying " );
		
		logger.debug( " writing patches " );
		ImgOps.writeFloat( ePatchImg, 
				String.format("/groups/jain/home/bogovicj/projects/crackPatching/edgelMatchPatch/patch_%03d.tif",i));
		ImgOps.writeFloat( fPatchImg, 
				String.format("/groups/jain/home/bogovicj/projects/crackPatching/edgelMatchPatch/patch_%03d_match.tif",i));
		logger.debug( " done writing patches " );
	}
	
	public void debugVisEdgelMatches(int i, Edgel e, List<Edgel> matches){
		
		String outFn = "/groups/jain/home/bogovicj/projects/crackPatching/edgelMatching";
		
		ArrayImgFactory<UnsignedByteType> ubfactory = new ArrayImgFactory<UnsignedByteType>();
		Img<UnsignedByteType> edgelMatchImg = ubfactory.create(mask, new UnsignedByteType());
		RandomAccess<UnsignedByteType> emiRa = edgelMatchImg.randomAccess();
		
		double[] pos = new double[e.numDimensions()];
		e.localize(pos);
		
		emiRa.setPosition( ArrayUtil.toIntRound(pos));
		emiRa.get().set( 255 );
		
		for ( Edgel match : matches )
		{
			match.localize(pos);
			emiRa.setPosition( ArrayUtil.toIntRound(pos));
			emiRa.get().set( 128 );
		}
		
//		logger.debug( "nnz: " + ImgOps.numNonZero(edgelMatchImg) );
		
		logger.debug( " writing match " );
		ImgOps.writeByte(edgelMatchImg, outFn + "/edgelMatches_" + i + ".tif");
		logger.debug( " done writing " );
	}
	
//	/**
//	 *
//	 */
//	public void computeAllAffinitiesIdxs()
//	{
//		
//		edgelTree = new KDTree<Edgel>( edgels, edgels );
//		switch ( search )
//		{
//		case COUNT:
//			countSearch = new  KNearestNeighborSearchOnKDTree<Edgel>( edgelTree, edgelSearchCount);
//			
//		default:
//			radiusSearch = new  RadiusNeighborSearchOnKDTree<Edgel>( edgelTree );
//
//		}
//		
//		edgelAffinities = new HashMap<EdgelIdxPair,Double>();
//
//		int i = 0;
//		for (Edgel e : edgels )
//		{
//
//			ArrayList<Integer> candidateEdgels = candidateEdgels(e);
//
//			for( Integer j : candidateEdgels )
//			{
//				EdgelIdxPair pair = new EdgelIdxPair(i,j);
//			
//				// if we already have this affinity skip the computation
//				if( edgelAffinities.containsKey(pair) ) continue;	
//
//				Edgel e2 = edgels.get(j);
//				double affinity = edgelAffinities( e, e2 ); 
//
//				edgelAffinities.put(pair, affinity);
//
//			}
//
//			i++;
//		}
//
//	}

	/**
	 * A pair of edgels
	 */ 
	public class EdgelPair
	{
		Edgel i;
		Edgel j;
		public EdgelPair(Edgel i, Edgel j)
		{
			if( compare(i,j) < 0 )
			{
				this.i = i;
				this.j = j;
			}
			else
			{
				this.i = j;
				this.j = i;
			}
		}
		/**
		 * provide a 
		 * @param i an edgel
		 * @param j another edgel
		 * @return -1 if i<j, 1 if i>j and zero if i==j
		 */
		public int compare(Edgel i, Edgel j){
			
			int nd = i.numDimensions();
			double[] ipos = new double[ i.numDimensions() ];
			i.localize(ipos);
			double[] jpos = new double[ j.numDimensions() ];
			j.localize(jpos);
			for (int d=0; d<nd; d++)
			{
				if( ipos[d] < jpos[d] ){
					return -1;
				}else if ( ipos[d] > jpos[d] ){
					return 1;
				}
				
			}
			return 0;
		}

	}
	
	/**
	 * A pair of indices
	 */ 
	public class EdgelIdxPair
	{
		int i;
		int j;
		public EdgelIdxPair(int i, int j)
		{
			if( i <= j)
			{
				this.i = i;
				this.j = j;
			}
			else
			{
				this.i = j;
				this.j = i;
			}
		}

		public int hashCode(){
	 		return j + i*numEdgels;	
		}
	}

}
