package net.imglib2.algorithms.crack;

import ij.IJ;

import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import net.imglib2.RandomAccess;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.algorithms.edge.EdgelTools;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;
import net.imglib2.util.LinAlgHelpers;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import edu.jhu.ece.iacl.utility.ArrayUtil;

public class EdgelClusteringRansac <T extends RealType<T> & NativeType<T>> 
{
	
	private static final float NEG = -1;
	private static final float POS =  1;
	
	private static final byte EMPTY    =  0;
	private static final byte TOUCHPOS =  1;
	private static final byte SEEDPOS  =  2;
	private static final byte TOUCHNEG = -1;
	private static final byte SEEDNEG  = -2;
	
	private final List<Edgel> edgels;
	private EdgelMatching<?> matcher;
	private int ndims;
	
	private float[] labelMems;
	private float   propMul 	= 1f;
	private float   eps     	= 0.000001f;
	
	private float   propThresh 	= 0.0001f;
	private float   dotThresh   = 0.5f;
	private   int   fitIters    = 5;
	
	private int maxIters = 3;
	
	private boolean success = false;
	
	LinkedList<Edgel> posEdgels;
	LinkedList<Edgel> negEdgels;
	
	ArrayList<Edgel> 	matches;
	HashSet <Edgel>		consSet;
	
	boolean[] 	seeded;
	int   		numNotSeed;
	
	private RandomDataGenerator rand = new RandomDataGenerator ();
	static Logger logger = LogManager.getLogger( EdgelClusteringRansac.class.getName() );
	

	public EdgelClusteringRansac(EdgelMatching<T> matcher){
		this.matcher = matcher;
		edgels = matcher.getEdgels();
		ndims = edgels.get(0).numDimensions();
		
		rand.reSeed( 316497258L );
	}
	
	public void setDotTolerance( float dotTol ){
		dotThresh = dotTol;
	}
	
	public void cluster()
	{
		
		int N = edgels.size();
		if( labelMems == null ){
			labelMems = new float[ N ];
		}
		logger.debug("clustering " + N + " edgels");
		
		numNotSeed = N;
		seeded = new boolean[ N ]; 
		
		posEdgels = new LinkedList<Edgel>();
		negEdgels = new LinkedList<Edgel>();
		consSet   = new HashSet<Edgel>();
		
		Edgel[] startPair = getStartingEdgels( 500 );
		Edgel fpos = startPair[0];
		Edgel fneg = startPair[1];
		
		// first pos
		labelMems[edgels.indexOf(fpos)] = 1f;
		double[] or = new double[ fpos.numDimensions() ];
		double[] dots = fitOrientationModel( or, fpos );
		boolean contradiction = updateConsensusSetMemsFirst( or, dots );
		
		// first neg
		labelMems[ edgels.indexOf(fneg) ] = -1f;
		Arrays.fill(or, 0);
		dots = fitOrientationModel( or, fneg );
		contradiction = updateConsensusSetMemsFirst( or, dots );
		
		logger.debug("pos edgel size " + posEdgels.size() );
		logger.debug("neg edgel size " + negEdgels.size() );
		
		logger.debug( " or " + ArrayUtil.printArray(or) ) ;
		logger.debug( " cset size " + consSet.size() ) ;
		logger.debug(" is contradiction: " + contradiction );
		
		
		runFirstIteration();
		
//		for ( int iter=0; iter<3; iter++)
//		{
//			logger.info("iteration: " + iter);
//			runIteration();
//		}
	}
	
	
	public void runFirstIteration()
	{
		double[] or   = new double[ ndims ];
		double[] dots = new double[ ndims ];
		
		boolean contradiction = false;
		
		int n = 0;
//		while(  posEdgels.size() > 0 && 
//				negEdgels.size() > 0 && 
//				n < 20000)
		int lastNotSeedCount = -1;
		int nscReps = 0;
		while( numNotSeed > 0 )
		{
			logger.debug( "nns " + numNotSeed );
			
			// do a positive edgel
			if( posEdgels.size() > 0 ){
				Edgel e = posEdgels.remove();
				Arrays.fill(or, 0);
				dots = fitOrientationModel( or, e );
				contradiction = updateConsensusSetMemsFirst( or, dots );
				logger.debug(" is contradiction: " + contradiction );
			}
			
			// do a negative edgel
			if( negEdgels.size() > 0 ){
				Edgel e = negEdgels.remove();
				Arrays.fill(or, 0);
				dots = fitOrientationModel( or, e );
				contradiction = updateConsensusSetMemsFirst( or, dots );
			}
			
			logger.trace(" is contradiction: " + contradiction );
			logger.trace(" pos edgel size " + posEdgels.size() );
			logger.trace(" neg edgel size " + negEdgels.size() );
			
			if(lastNotSeedCount == numNotSeed  ){
				nscReps++;
			}else{
				nscReps = 0;
			}
			
			if ( nscReps == 10 ){
				break;
			}
			
			lastNotSeedCount = numNotSeed;
			
			n++;
		}
	}
	
	public Edgel[] getStartingEdgels ( int maxTries ) 
	{
		Edgel[] firstMatches = new Edgel[2];
		
		int n = 0;
		while( n < maxTries ){
			
			int i = rand.nextInt(0, edgels.size()-1); 
			Edgel e = edgels.get(i);
			
			Edgel f = findOpposite( e );
			
			if( f != null ){
				firstMatches[0] = e;
				firstMatches[1] = f;
				break;
			}
			n++;
		}
		logger.debug("Found starting edgels after " + n + " iterations.");
		return firstMatches;
	}
	
	public Edgel findOppositeLoopy( Edgel e ){
		int best = -1;
		boolean first = true;
		while( best < 0 && matcher.getEdgelSearchCount() < 61 )
		{
			if( first ){ first = false; }
			else{ 
				matcher.setEdgelSearchCount( 
					matcher.getEdgelSearchCount() + 2 ); 
			}
			
			ArrayList<Edgel> matches = matcher.candidateEdgels(e);
			double[] dots = matcher.dotProduct(e, matches);
			double min = -0.01f;
			for( int i = 0; i<dots.length; i++){
				if( dots[i] < min ){
					min = dots[i];
					best = i;
				}
			}
		}
		return matches.get(best);
	}
	
	public Edgel findOpposite( Edgel e ){
		int best = -1;

		ArrayList<Edgel> matches = matcher.candidateEdgels(e);
		double[] dots = matcher.dotProduct(e, matches);
		double min = -0.01f;
		for( int i = 0; i<dots.length; i++){
			if( dots[i] < min ){
				min = dots[i];
				best = i;
			}
		}
		
		if ( best == -1 ){ return null; }
		
		return matches.get(best);
	}

	public double[] fitOrientationModel( double[] orientation, Edgel e ){
		
		int i = edgels.indexOf(e);
		logger.trace(" fitting model around edgel idx " + i);
		
		if(  !seeded[ i ] )
		{
			seeded[ i ] = true;
			numNotSeed--;
		}
		
		matches = matcher.candidateEdgels(e);
		matches.add(e);
		
		consSet.clear();
		consSet.addAll( matches );
		consSet.add(e);
		
		double[] dots = null;
		
		for (int iter=0; iter < fitIters; iter++)
		{
			boolean anyChanges = false;
			
			logger.trace(" iter: " + iter );
			
			// estimate orientation from consensus set
			EdgelTools.averageGradientDirectionNaive( consSet, orientation );

			logger.trace( " or " + ArrayUtil.printArray(orientation) );
			
			// update consensus set
			int N = matches.size();
			dots = matcher.dotProduct( orientation, matches );
			
			for( i=0; i<N; i++ ) 
			{	
				if( dots[i] > dotThresh ) { 
					boolean wasChange = consSet.add( matches.get(i) );
					anyChanges = anyChanges || wasChange;
				}else{
					boolean wasChange = consSet.remove(  matches.get(i) );
					anyChanges = anyChanges || wasChange;
				}
				
			}
			
			// exit if consensus set did not change
			if( !anyChanges ){
				logger.trace(" exit after " + iter + " iters");
				break;
			}
		}
		
		return dots;
	}
	
	public boolean updateConsensusSetMemsFirst( double[] or, double[] dots )
	{
		logger.trace("Updating memberships, cset of size: " + consSet.size() );
		
		boolean first = true;
		boolean contradiction = false;
		
		float val = Float.NaN;
		double posSum = 0;
		double negSum = 0;
		
		int[] idxs = new int[ consSet.size() ];
		
		int i = 0;
		for( Edgel c : consSet )
		{	
			int j = edgels.indexOf( c );
			idxs[i] = j;
			
			if(  first && 
				( labelMems[j] >  propThresh || 
				  labelMems[j] < -propThresh ))
			{
				val = labelMems[ j ];
				first = false;
			}
			else if( labelMems[j] * val < 0 )
			{
				contradiction = true;
			}
			
			if( labelMems[j] > propThresh ){
				posSum += labelMems[j];
			}else if(  labelMems[j] < -propThresh ){
				negSum -= labelMems[j]; // make negsum positive
			}
			
			i++;
		}
		
		double tmp = posSum + negSum;
		if( tmp > eps )
		{
			posSum = posSum / tmp;
			negSum = negSum / tmp;
		}
		
		logger.trace("posSum " + posSum );
		logger.trace("negSum " + negSum );
		
		if( posSum > negSum ){
			val =  (float)posSum;
		}else{
			val = -(float)negSum;
		}
		
		for( i=0; i<consSet.size(); i++ )
		{
			int j = idxs[i];
			labelMems[j] = (float)( 0.7 * val + 0.3 * labelMems[j]); 
			
			logger.trace("labelMem " + j + " " + labelMems[j] );
			
			if( !seeded[j] )
			{
				if( labelMems[j] > 0  && 
					!posEdgels.contains( edgels.get(j)))
				{
					posEdgels.add( edgels.get(j) );
				}
				else if( labelMems[j] < 0 && 
						 !negEdgels.contains( edgels.get(j)))
				{
					negEdgels.add( edgels.get(j) );
				}
			}
		}
		
		return contradiction;
	}
	
	public boolean isContradiction( double[] or )
	{
		boolean first = true;
		float val = Float.NaN;
		
		for( Edgel c : consSet )
		{	
			int j = edgels.indexOf( c );
			
			if( first && labelMems[j] != 0 )
			{
				val = labelMems[ j ];
				first = false;
			}
			else if( labelMems[j] * val < 0 )
			{
				return true;
			}
		}
		return false;
	}

	public void makeEdgelClusterImgMem( Img<T> clusterImg ){
		
		RandomAccess<T> ra = clusterImg.randomAccess();
		
		int N = edgels.size();
		for ( int i = 0; i < N; i++ )
		{
			Edgel e = edgels.get(i);
			for( int d=0; d<clusterImg.numDimensions(); d++){
				ra.setPosition( (int)(Math.round( e.getDoublePosition(d) )), d);
			}
			ra.get().setReal(labelMems[i]);
			
		}
	}

	private int nextUnlabeled( int start ){
		for ( int i=start; i<edgels.size(); i++){
			if( labelMems[i] == 0 ) return i;
		}
		return -1;
	}
	private int nextUnlabeled(  ){
		RandomDataGenerator rand = new RandomDataGenerator();
		int start = rand.nextInt(0, edgels.size()-1);
		for ( int i=start; i<edgels.size(); i++){
			if( labelMems[i] == 0 ) return i;
		}
		return -1;
	}

	
	public static void main(String[] args) {
		
		int downSampleFactor = 4;
//		double searchRadius = 100;
		int    searchCount  = 60;
		
		String imgfn = "/data-ssd1/john/projects/crackSegmentation/groundTruth/closeup/img_ds"+downSampleFactor+".tif";
		String maskfn = "/data-ssd1/john/projects/crackSegmentation/groundTruth/closeup/labels_interp_smooth_ds"+downSampleFactor+".tif";
		Img<FloatType> img =  ImagePlusAdapter.convertFloat( IJ.openImage(imgfn) );
		Img<FloatType> mask = ImagePlusAdapter.convertFloat( IJ.openImage(maskfn) );
		
//		int[] patchSize = new int[] { 45, 45, 19 };
		int[] patchSize = new int[] { 31, 31, 19 };
//		int[] patchSize = new int[] { 15, 15, 13 };
		
		CrackCorrection<FloatType> cc = new CrackCorrection<FloatType>(
				img, mask, patchSize);
		
		EdgelMatching<FloatType> em = new EdgelMatching<FloatType>(img, mask, patchSize);
		em.debugDir = "/data-ssd1/john/projects/crackPatching/cropEdgelMatch";
		em.debugSuffix = "ds"+downSampleFactor;
		
		em.setSearchType( EdgelMatching.SearchTypes.COUNT );
		em.setEdgelSearchCount( searchCount );
		
//		em.setEdgelSearchRadius( searchRadius / (downSampleFactor) );
		
		cc.edgelMatcher = em;
		
		cc.computeEdgels();
		cc.edgelMatcher.setEdgels(cc.getEdgels());
		
		EdgelClusteringRansac<FloatType> ec = new EdgelClusteringRansac<FloatType>(
			cc.edgelMatcher);
		
		ec.cluster( );
		
		Img<FloatType> clusterImg = img.factory().create( img, img.firstElement());
//		ec.makeEdgelClusterImg( clusterImg );
		ec.makeEdgelClusterImgMem( clusterImg );
		
		Format formatter = new SimpleDateFormat("yyyyMMdd-HHmmss");
		String dtstr = formatter.format(Calendar.getInstance().getTime());
		String fn = String.format("/groups/saalfeld/home/bogovicj/tmp/edgelClustersMemRansac_%s.tif", dtstr);
		logger.info("fn: " + fn);
		ImgOps.writeFloat( clusterImg, fn);
		
		
		logger.info("finished");
		System.exit(0);
	}

}
