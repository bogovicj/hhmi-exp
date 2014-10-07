package net.imglib2.algorithms.opt.astar;

import java.util.ArrayList;

import net.imglib2.algorithms.patch.SubPatch2dLocation;


public class AStarMax<T extends Comparable<T>> extends AbstractAStar<T> {
	
	public AStarMax(){
		super();
	}
	
	public AStarMax( T rootData ){
		super(rootData);
	}
	
//	@Override
//	protected boolean isFinished() {
//		// TODO Auto-generated method stub
//		return false;
//	}
//
//	@Override
//	protected boolean continueDepthSearch() {
//		// TODO Auto-generated method stub
//		return true;
//	}
	
	@Override
	protected boolean propagateFunction(SortedTreeNode<T> currentBest, SortedTreeNode<T> child) {
		boolean shouldParentChange = true;
		if( child.getData().compareTo( currentBest.getData()) < 0 ){
			child.setData( currentBest.getData() );
			shouldParentChange = false;
		}
		return shouldParentChange;
	}
	
	public static void main( String[] args ){
		logger.info("Starting");
		
		SubPatch2dLocation  root = new SubPatch2dLocation(0, 0, 0.0);
		AStarMax<SubPatch2dLocation> astar = new AStarMax<SubPatch2dLocation>(root);
		
		ArrayList<SubPatch2dLocation> layer1 = new ArrayList<SubPatch2dLocation>();
		layer1.add( new SubPatch2dLocation( 1, 1, 4.0 ));
		layer1.add( new SubPatch2dLocation( 2, 1, 2.0 ));
		layer1.add( new SubPatch2dLocation( 1, 3, 6.0 ));
		layer1.add( new SubPatch2dLocation( 2, 3, 8.0 ));
		
		astar.addCandidates( layer1 );
		
		logger.info("best 1: " + astar.getBest() );
		
		ArrayList<SubPatch2dLocation> layer2 = new ArrayList<SubPatch2dLocation>();
		layer2.add( new SubPatch2dLocation( 1, 1,  6.0 ));
		layer2.add( new SubPatch2dLocation( 2, 1,  1.0 ));
		layer2.add( new SubPatch2dLocation( 1, 3,  6.0 ));
		layer2.add( new SubPatch2dLocation( 2, 3, 10.0 ));
		
		astar.addCandidates( layer2 );
		
		logger.info("best 2: " + astar.getBest() );
		
		/*
		SortedTreeNode<Double> root = new SortedTreeNode<Double>(new Double(5));
		logger.info("root: " + root );
		
		root.addChild( new SortedTreeNode<Double>(new Double( 6 )) );
		root.addChild( new SortedTreeNode<Double>(new Double( 4 )) );
		
		logger.info("root: " + root );
		
		Set<SortedTreeNode<Double>> children = root.getChildren();
		
		Iterator<SortedTreeNode<Double>> it = children.iterator();
		logger.info("child 1: " + it.next());
		logger.info("child 2: " + it.next());
		*/
		
		logger.info("Finished");
	}




}
