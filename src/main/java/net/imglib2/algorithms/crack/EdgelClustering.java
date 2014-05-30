package net.imglib2.algorithms.crack;

import ij.IJ;

import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.List;

import net.imglib2.RandomAccess;
import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgOps;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

import edu.jhu.ece.iacl.utility.ArrayUtil;

public class EdgelClustering <T extends RealType<T> & NativeType<T>> 
{
		
	public static final byte UNKNOWN   =  0;
	public static final byte CONFLICT  =  3;
	
	public static final byte CLASSA =  1;
	public static final byte CLASSB =  2;
	
	
	private final List<Edgel> edgels;
	private EdgelMatching<?> matcher;
	
	private boolean[][] nbrTable;
	private byte[] labels;
	
	
	private float[] labelMems;
	private float   propMul 	= 1f;
	private float   eps     	= 0.000001f;
	
	private float   propThresh 	= 0.1f;
	private float   dotThresh   = 0.3f;
	
	private int maxIters = 3;
	
	static Logger logger = LogManager.getLogger( EdgelClustering.class.getName() );
	
	public EdgelClustering( List<Edgel> edgels, EdgelMatching<T> matcher){
		this.edgels = edgels;
		this.matcher = matcher;
	}
	
	public void cluster(){
		
//		findNeighbors();
		
		int nConflicts        = 9999;
		int lastIterConflicts = 9999;
		
		for( int i = 0; i<maxIters; i++ )
		{

			logger.info("iteration " + i + " of " + maxIters );
			
			// run an iteration
//			nConflicts = iteration();
			nConflicts = iterationMem();
			
//			if ( nConflicts == 0 ){
//				break;
//			}
//			
//			if ( lastIterConflicts == nConflicts ){
//				break;
//			}
			
			lastIterConflicts = nConflicts;
			
			ArrayUtil.multiply(labelMems, 2);
		}
	}
	
	public void makeEdgelClusterImg( Img<T> clusterImg ){
		
		RandomAccess<T> ra = clusterImg.randomAccess();
		
		int N = edgels.size();
		for ( int i = 0; i < N; i++ )
		{
			Edgel e = edgels.get(i);
			for( int d=0; d<clusterImg.numDimensions(); d++){
				ra.setPosition( (int)(Math.round( e.getDoublePosition(d) )), d);
			}
			if( labels[i] == CLASSA ){
				ra.get().setReal( 1 );
			}else if( labels[i] == CLASSB ){
				ra.get().setReal( 2 );
			}else if( labels[i] == CONFLICT ){
				ra.get().setReal( 3 );
			}
			
		}
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
	
	public int iterationMem()
	{
		int N = edgels.size();
		if( labelMems == null ){
			labelMems = new float[ N ];
		}
		
		for ( int i = 0; i < N; i++ )
//		for ( int i = 0; i < 3; i++ )
		{
			if ( i % 5000 == 0){
				logger.info("   edgel " + i + " of " + edgels.size() );
			}
			
			Edgel e = edgels.get(i);
			ArrayList<Edgel> matches = matcher.candidateEdgels( e );
//			boolean[] sameSide = matcher.similarlyOriented(e, matches);
			float[] dots = matcher.dotProduct(e, matches);

			// initialize
			if( i == 0){ labelMems[i] = 1.0f; }
			
			float thisClass = labelMems[i];
			
			if( Math.abs(thisClass) < eps ){
//				logger.info(" find best among " + matches.size() + " matches ");
//				labelMems[i] = bestClassMem( sameSide, matches );
				labelMems[i] = bestClassMem( dots, matches );
			}else if ( Math.abs(thisClass) > propThresh ){
//				logger.info(" propagate to " + matches.size() + " matches ");
//				propagateClassMem( i, sameSide, matches );
				propagateClassMem( i, dots, matches );
			}
			
		}
		
		return 1;
	}
	public float bestClassMem( float[] dots, List<Edgel> matches ){

		float memOut = 0.0f;
		float dotSum = 0.0f;
		for ( int j=0; j<matches.size(); j++ )
		{
			Edgel match = matches.get(j);
			int k = edgels.indexOf(match);

			if ( Math.abs(labelMems[k]) > propThresh &&
				 Math.abs(dots[j]) > dotThresh	)
			{
				memOut += dots[j];
			}
		}
		memOut /= matches.size();

		return memOut;
	}
	
	public void propagateClassMem( int i, float[] dots, List<Edgel> matches ){
		
		float thisMem = labelMems[i];
		
		for ( int j=0; j<matches.size(); j++ )
		{
			Edgel match = matches.get(j);
			int k = edgels.indexOf(match);
			
			if( Math.abs(labelMems[k]) < eps &&
				Math.abs(dots[j]) > dotThresh )
			{
				labelMems[k] =   propMul * dots[j] * thisMem;
			}
		}
	}

	
	public void propagateClassMem( int i, boolean[] sameSide, List<Edgel> matches ){
		
		float thisMem = labelMems[i];
		
		for ( int j=0; j<matches.size(); j++ )
		{
			Edgel match = matches.get(j);
			int k = edgels.indexOf(match);
			
			if( labelMems[k] < eps ){
				if( sameSide[j] ){
					labelMems[k] =   propMul * thisMem;
				}else{
					labelMems[k] = - propMul * thisMem;
				}
			}
		}
	}
	
	public int iteration()
	{
		int N = edgels.size();
		if( labels == null ){
			labels = new byte[ N ];
		}
		Arrays.fill(labels, UNKNOWN);
		
		for ( int i = 0; i < N; i++ )
//		for ( int i = 0; i < 3; i++ )
		{
			if ( i % 5000 == 0){
				logger.info("   edgel " + i + " of " + edgels.size() );
			}
			
			Edgel e = edgels.get(i);
			ArrayList<Edgel> matches = matcher.candidateEdgels( e );
			boolean[] sameSide = matcher.similarlyOriented(e, matches);

			// initialize
			if( i == 0){ labels[i] = CLASSA; }
			
			byte thisClass = labels[i];
			
			if( thisClass == UNKNOWN ){
//				logger.info(" find best among " + matches.size() + " matches ");
				labels[i] = bestClass( sameSide, matches );
			}else{
//				logger.info(" propagate to " + matches.size() + " matches ");
				propagateClass( i, sameSide, matches );
			}
			
		}
		
		return 1;
	}

	
	public void propagateClass( int i, boolean[] sameSide, List<Edgel> matches ){
		
		byte thisClass = labels[i];
		byte othrClass = otherClass( thisClass );
		
		for ( int j=0; j<matches.size(); j++ )
		{
			Edgel match = matches.get(j);
			int k = edgels.indexOf(match);
			
			if( labels[k] == UNKNOWN ){
				if( sameSide[j] ){
					labels[k] = thisClass;
				}else{
					labels[k] = othrClass;
				}
			}
//			else{
//				if( isConflict(  thisClass, labels[k], sameSide[j]) ){
//					labels[k] = CONFLICT;
//					labels[i] = CONFLICT;
//				}
//			}
		}
	}
	
	
	public byte otherClass( byte thisClass ){
		if( thisClass == CLASSA ){
			 return CLASSB;
		}else if( thisClass == CLASSB){
			return CLASSA;
		}
		return UNKNOWN;
	}
	
	
	
	public byte bestClass( boolean[] sameSide, List<Edgel> matches ){
		
		int aCount = 0;
		int bCount = 0;
		
		for ( int j=0; j<matches.size(); j++ )
		{
			Edgel match = matches.get(j);
			int k = edgels.indexOf(match);
			
			if( labels[k] == CLASSA ){
				if( sameSide[j] ){
					aCount++;	
				}else{
					bCount++;
				}
				
			}else if( labels[k] == CLASSB ){
				if( sameSide[j] ){
					bCount++;
				}else{
					aCount++;
				}
			}
		}
		
//		if( aCount > 0 && bCount == 0 ){
//			return CLASSA;
//		}else if( bCount > 0 && aCount == 0){
//			return CLASSB;
//		}
		if( aCount > 2*bCount ){
			return CLASSA;
		}else if( bCount > 2*aCount ){
			return CLASSB;
		}
		return UNKNOWN;
	}
	
	/**
	 * Returns true if htere is a conflict between the assigned labels and
	 * the sameSide flag.
	 * @param label1
	 * @param label2
	 * @param sameSide
	 * @return
	 */
	private static boolean isConflict( byte label1, byte label2, boolean sameSide){
		return ( ( label1 == label2 ) ^ sameSide );
	}
	
	public void findNeighbors()
	{
		int N = edgels.size();
		
		nbrTable = new boolean[N][N];
		ArrayUtil.fill(nbrTable, false);
		
		int i = 0;
		for( Edgel e : edgels )
		{	
			ArrayList<Edgel> matches = matcher.candidateEdgels( e );
			boolean[] sameSide = matcher.similarlyOriented(e, matches);
			
			if ( i % 200 == 0){
				logger.info(" edgel " + i + " of " + edgels.size() );
			}
			for ( int j=0; j<matches.size(); j++ )
			{
				Edgel match = matches.get(j);
				int k = edgels.indexOf(match);
				
				nbrTable[i][k] = sameSide[j];
			}
			
			i++;
		}
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
		
		EdgelClustering<FloatType> ec = new EdgelClustering<FloatType>(
				cc.getEdgels(), cc.edgelMatcher);
		
		ec.cluster( );
		Img<FloatType> clusterImg = img.factory().create( img, img.firstElement());
//		ec.makeEdgelClusterImg( clusterImg );
		ec.makeEdgelClusterImgMem( clusterImg );
		
		Format formatter = new SimpleDateFormat("yyyyMMdd-HHmmss");
		String dtstr = formatter.format(Calendar.getInstance().getTime());
		String fn = String.format("/groups/jain/home/bogovicj/tmp/edgelClustersMem_%s.tif", dtstr);
		logger.info("fn: " + fn);
		ImgOps.writeFloat( clusterImg, fn);
		
		
		/* CHECK CONFLICT CONDITION METHOD */
//		System.out.println( "A A S :"  + isConflict( CLASSA, CLASSA, true));
//		System.out.println( "A B S :"  + isConflict( CLASSA, CLASSB, true));
//		System.out.println( "A A D :"  + isConflict( CLASSA, CLASSA, false));
//		System.out.println( "A B D :"  + isConflict( CLASSA, CLASSB, false));
		
		
		logger.info("finished");
		System.exit(0);
	}

}
