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
	
	private float[] labelMems;
	private float   propMul 	= 1f;
	private float   eps     	= 0.000001f;
	
	private float   propThresh 	= 0.0001f;
	private float   dotThresh   = 0.5f;
	private   int   fitIters    = 7;
	
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
	}
	
	public void setDotTolerance( float dotTol ){
		dotThresh = dotTol;
	}
	
	public void cluster(){
		
//		findNeighbors();
		
		int N = edgels.size();
		if( labelMems == null ){
			labelMems = new float[ N ];
		}
		
		numNotSeed = N;
		seeded = new boolean[ N ]; 
		
		posEdgels = new LinkedList<Edgel>();
		negEdgels = new LinkedList<Edgel>();
		
		Edgel match = findOpposite( edgels.get(0) );

		logger.debug( " " + edgels.get(0) );
		logger.debug( " " + match ) ;
		
		// first pos
		labelMems[0] = 1f;
		double[] or = new double[ edgels.get(0).numDimensions() ];
		double[] dots = fitOrientationModel( or, edgels.get(0) );
		boolean contradiction = updateConsensusSetMems( or, dots );
		
		// first neg
		labelMems[ edgels.indexOf(match) ] = -1f;
		Arrays.fill(or, 0);
		dots = fitOrientationModel( or, match );
		contradiction = updateConsensusSetMems( or, dots );
		
		System.out.println(" pos edgel size " + posEdgels.size() );
		System.out.println(" neg edgel size " + negEdgels.size() );
		
		logger.debug( " or " + ArrayUtil.printArray(or) ) ;
		logger.debug( " cset size " + consSet.size() ) ;
		logger.debug(" is contradiction: " + contradiction );
		
		
//		while( numNotSeed > 0 )
//		{
//			// do a positive edgel
//			Edgel e = posEdgels.remove();
//			Arrays.fill(or, 0);
//			dots = fitOrientationModel( or, e );
//			contradiction = updateConsensusSetMems( or, dots );
//			
//			// do a negative edgel
//			e = negEdgels.remove();
//			Arrays.fill(or, 0);
//			dots = fitOrientationModel( or, e );
//			contradiction = updateConsensusSetMems( or, dots );
//			
//			System.out.println( "nns " + numNotSeed );
//		}
		
		
		
//		for( int i = 0; i<maxIters; i++ )
//		{
//
//			logger.info("iteration " + i + " of " + maxIters );
//			
//			// run an iteration
//			
//			lastIterConflicts = nConflicts;
//			
//		}
	}
	
	public Edgel findOpposite( Edgel e ){
		ArrayList<Edgel> matches = matcher.candidateEdgels(e);
		double[] dots = matcher.dotProduct(e, matches);
		double min = -0.1f;
		int best = -1;
		for( int i = 0; i<dots.length; i++){
			if( dots[i] < min ){
				min = dots[i];
				best = i;
			}
		}
		return matches.get(best);
	}
	
	public void firstIter(){
		
	}
	
	public double[] fitOrientationModel( double[] orientation, Edgel e ){
		
		seeded[ edgels.indexOf(e) ] = true;
		numNotSeed--;
		
		matches = matcher.candidateEdgels(e);
		
		consSet = new HashSet<Edgel>();
		consSet.addAll( matches );
		
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
			
			
			for( int i=0; i<N; i++ ) 
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
			logger.trace(" consensus set size: " + consSet.size());
			
			if( !anyChanges ){
				logger.trace(" exit after " + iter + " iters");
				break;
			}
		}
		
		return dots;
	}
	
	public void updateConsensusSetMemsOld( double[] dots, double[] or )
	{		
		double max = 0;
		
		for( Edgel c : consSet )
		{	
			int j = edgels.indexOf( c );
			
			if( labelMems[j] != 0 ){
				labelMems[j] = (float)( LinAlgHelpers.dot( or, c.getGradient()) );
			}
		}
	}
	
	public boolean updateConsensusSetMems( double[] or, double[] dots )
	{
		logger.debug("Updating memberships, cset of size: " + consSet.size() );
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
			
			if( first && labelMems[j] != 0 )
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
				negSum -= labelMems[j];
			}
			
			i++;
		}
		
		double tmp = posSum + negSum;
		if( tmp > eps ){
			posSum = posSum / tmp;
			negSum = negSum / tmp;
		}
		
		if( posSum > negSum ){
			val =  (float)posSum;
		}else{
			val = -(float)negSum;
		}
		
		for( i=0; i<consSet.size(); i++ )
		{
			int j = idxs[i];
			labelMems[j] = (float)( 0.7 * val + 0.3 * labelMems[j]); 
			
			if( labelMems[j] > 0 && !seeded[j])
			{
				posEdgels.add( edgels.get(j) );
			}
			else if( labelMems[j] > 0 && !seeded[j] )
			{
				negEdgels.add( edgels.get(j) );
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
