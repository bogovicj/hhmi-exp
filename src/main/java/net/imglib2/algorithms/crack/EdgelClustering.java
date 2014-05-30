package net.imglib2.algorithms.crack;

import java.util.ArrayList;
import java.util.List;

import net.imglib2.algorithm.edge.Edgel;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

public class EdgelClustering <T extends RealType<T> & NativeType<T>> {
	
	static Logger logger = LogManager.getLogger( EdgelClustering.class.getName() );
	
	private final List<Edgel> edgels;
	private EdgelMatching<?> matcher;
	
	private int[][] classTable;
	
	public EdgelClustering( List<Edgel> edgels, EdgelMatching<T> matcher){
		this.edgels = edgels;
		this.matcher = matcher;
	}
	
	public void cluster(){
		int N = edgels.size();
		classTable = new int[N][N];
		int i = 0;
		for( Edgel e : edgels ){
			
			ArrayList<Edgel> matches = matcher.candidateEdgels( e );
			matcher.filterEdgelsByNormal(e, matches);
			
			i++;
		}
		
	}
	
	public static void main(String[] args) {
		
		logger.info("finished");
		System.exit(0);
	}

}
