package net.imglib2.algorithms.crack;

import ij.IJ;

import java.util.ArrayList;
import java.util.Arrays;
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
	public static final byte CONFLICT  = -2;
	
	public static final byte CLASSA = -1;
	public static final byte CLASSB =  1;
	
	
	private final List<Edgel> edgels;
	private EdgelMatching<?> matcher;
	
	private boolean[][] nbrTable;
	private byte[] labels;
	
	private int maxIters = 1;
	
	static Logger logger = LogManager.getLogger( EdgelClustering.class.getName() );
	
	public EdgelClustering( List<Edgel> edgels, EdgelMatching<T> matcher){
		this.edgels = edgels;
		this.matcher = matcher;
	}
	
	public void cluster(){
		
		findNeighbors();
		
		int nConflicts        = 9999;
		int lastIterConflicts = 9999;
		
		for( int i = 0; i<maxIters; i++ )
		{
			
			// run an iteration
			nConflicts = iteration();
			
			if ( nConflicts == 0 ){
				break;
			}
			
			if ( lastIterConflicts == nConflicts ){
				break;
			}
			
			lastIterConflicts = nConflicts;
			
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
			ra.get().setOne();
			if( labels[i] == CLASSB ){
				ra.get().mul(2.0);
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
		{
			
			Edgel e = edgels.get(i);
			ArrayList<Edgel> matches = matcher.candidateEdgels( e );
			boolean[] sameSide = matcher.similarlyOriented(e, matches);

			// initialize
			if( i == 0){ labels[i] = CLASSA; }
			
			byte thisClass = labels[i];
			
			if( thisClass == UNKNOWN ){
				labels[i] = bestClass( sameSide, matches );
			}else{
				propagateClass( i, sameSide, matches );
			}
			
		}
		
		return 0;
	}
	
	public void propagateClass( int i, boolean[] sameSide, List<Edgel> matches ){
		
		byte thisClass = labels[i];
		byte othrClass = (byte)(-1 * labels[i]);
		
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
			}else{
				if( isConflict(  thisClass, labels[k], sameSide[j]) ){
					labels[k] = CONFLICT;
					labels[i] = CONFLICT;
				}
			}
		}
	}
	
	public byte bestClass( boolean[] sameSide, List<Edgel> matches ){
		int aCount = 0;
		int bCount = 0;
		
		for ( int j=0; j<matches.size(); j++ )
		{
			Edgel match = matches.get(j);
			int k = edgels.indexOf(match);
			
			if( labels[k] == CLASSA ){
				aCount++;
			}else if( labels[k] == CLASSB ){
				bCount++;
			}
		}
		
		if( aCount > 0 & bCount ==0 ){
			return CLASSA;
		}else if( bCount > 0 & aCount == 0){
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
			
			if ( i % 50 == 0){
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
		int    searchCount  = 50;
		
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
		ec.makeEdgelClusterImg( clusterImg );
		
		ImgOps.writeFloat( clusterImg, "/groups/jain/home/bogovicj/tmp/edgelClusters.tif");
		
//		System.out.println( "A A S :"  + isConflict( CLASSA, CLASSA, true));
//		System.out.println( "A B S :"  + isConflict( CLASSA, CLASSB, true));
//		System.out.println( "A A D :"  + isConflict( CLASSA, CLASSA, false));
//		System.out.println( "A B D :"  + isConflict( CLASSA, CLASSB, false));
		
		
		logger.info("finished");
		System.exit(0);
	}

}
