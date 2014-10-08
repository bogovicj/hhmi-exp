package net.imglib2.algorithms.patch;

import java.util.*;

import net.imglib2.algorithms.opt.astar.SortedTreeNode;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;

public class Patch2dFill3d< T extends Comparable<T> > {

	// logging
	public static Logger logger = LogManager.getLogger(
									Patch2dFill3d.class.getName());
	
	private TreeSet<SortedTreeNode<T>> set;
	private SortedTreeNode<T> root;
	
	public Patch2dFill3d( T rootDat ){
		root = new SortedTreeNode<T>(rootDat );
		set = new TreeSet<SortedTreeNode<T>>();
		set.add( root );
	}
	
	public SortedTreeNode<T> next(){
		return set.first();
	}
	
	public SortedTreeNode<T> getRoot(){
		return root;
	}
	
	public void addFromNode( SortedTreeNode<T> parent, List<T> newPatches ){
		addFromNode( parent, newPatches, true );
	}
	
	public void addFromNode( SortedTreeNode<T> parent, List<T> newPatches, boolean removeParents ){
		Iterator<T> it = newPatches.iterator();
		
		while( it.hasNext() ){
			T nextData = it.next();
			SortedTreeNode<T> child = new SortedTreeNode<T>( nextData, parent );
			parent.addChild(child);
			set.add( child );
		}
		
		if( removeParents )
			set.remove( parent );
	}
	
	
	public static void main(String[] args) {
		logger.info("Start");
		
		SubPatch2dLocation  root = new SubPatch2dLocation(0, 0, 0, 0.0);
		Patch2dFill3d<SubPatch2dLocation> pf = new Patch2dFill3d<SubPatch2dLocation>(root);
		
		ArrayList<SubPatch2dLocation> layer1 = new ArrayList<SubPatch2dLocation>();
		layer1.add( new SubPatch2dLocation( 1, 1, 1, 1.0 ));
		layer1.add( new SubPatch2dLocation( 2, 1, 2, 3.0 ));
		layer1.add( new SubPatch2dLocation( 1, 3, 3, 5.0 ));
		layer1.add( new SubPatch2dLocation( 2, 3, 4, 7.0 ));
		
		SortedTreeNode<SubPatch2dLocation> rootNode = pf.next();
		pf.addFromNode( rootNode, layer1 );
		
		SortedTreeNode<SubPatch2dLocation> nextBest = pf.next();
		logger.info("\nnext Best: " + nextBest );
		logger.info("next Best Parent: " + nextBest.getParent() );
		
		ArrayList<SubPatch2dLocation> layer2 = new ArrayList<SubPatch2dLocation>();
		layer2.add( new SubPatch2dLocation( 1, 1, 1, 2.0 ));
		layer2.add( new SubPatch2dLocation( 2, 1, 2, 4.0 ));
		layer2.add( new SubPatch2dLocation( 1, 3, 3, 6.0 ));
		layer2.add( new SubPatch2dLocation( 2, 3, 4, 7.0 ));
		
		pf.addFromNode( nextBest, layer2 );
		
		nextBest = pf.next();
		logger.info("\nnext Best: " + nextBest );
		logger.info("next Best Parent: " + nextBest.getParent() );
		logger.info("next Best Parent^2: " + nextBest.getParent().getParent() );
		
		logger.info("Finish");
	}

}