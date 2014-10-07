package net.imglib2.algorithms.opt.astar;


import java.util.*;

import org.apache.log4j.*;

public abstract class AbstractAStar <T extends Comparable<T>> {

	
	SortedTreeNode<T> root;
	SortedTreeNode<T> bestCandidate;
	
	// logging
	public static Logger logger = LogManager.getLogger(
									AbstractAStar.class.getName());
	
	public AbstractAStar(){
	}
	
	public AbstractAStar( T rootData ){
		setRoot( rootData );
	}
	
	public void setRoot( T rootData ){
		root = new SortedTreeNode<T>( rootData );
		bestCandidate = root;
	}
	
	public SortedTreeNode<T> getBest(){
		return bestCandidate;
	}
	
	
	public void addCandidates( Collection<T> nextLayer ){
		Iterator<T> it = nextLayer.iterator();
		while( it.hasNext() ){
			bestCandidate.addChild( 
					new SortedTreeNode<T>( it.next() ));
		}
		
		logger.info("best (in): " + bestCandidate );
		
		// the new best candidate is the best among the added candidates 
		logger.info("children: " + bestCandidate.getChildren() );
		
		SortedTreeNode<T> newBest = bestCandidate.getChildren().iterator().next();
		
		propagateUp( bestCandidate, newBest );
		
		bestCandidate = newBest;
		
	}
	
	
//	/**
//	 * Returns whether a local optimum has been found
//	 * @return
//	 */
//	protected abstract boolean isFinished();
//	
//	protected abstract boolean continueDepthSearch();
	
	protected abstract boolean propagateFunction( SortedTreeNode<T> currentBest, SortedTreeNode<T> child );
	
	protected void propagateUp( SortedTreeNode<T> currentBest, SortedTreeNode<T> child ){
		propagateFunction( currentBest, child );
		
	}
}