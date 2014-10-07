package net.imglib2.algorithms.opt.astar;

import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;


public class SortedTreeNode<T extends Comparable<T>> implements Comparable<SortedTreeNode<T>> {
	private T							data;
	private SortedTreeNode<T> 			parent;
	private SortedSet<SortedTreeNode<T>> children;
	
	/**
	 * Constructor
	 * @param data
	 * @param parent
	 */
	public SortedTreeNode(T data, SortedTreeNode<T> parent ){
		this.data = data;
		this.parent = parent;
		children = new TreeSet<SortedTreeNode<T>>();
	}
	
	/**
	 * Constructor for the root.
	 * @param data
	 */
	public SortedTreeNode( T data ){
		this( data, null );
	}
	
	public void addChild( SortedTreeNode<T> child ){
		child.setParent( this );
		children.add( child );
	}
	
	public void addChildren( Collection<SortedTreeNode<T>> children ){
		Iterator<SortedTreeNode<T>> it = children.iterator();
		while (it.hasNext()){
			addChild( it.next() );
		}
	}
	
	public T getData(){
		return data;
	}
	
	public void setData( T data){
		this.data = data;
	}
	
	protected void setParent( SortedTreeNode<T> parent ){
		this.parent = parent;
	}
	
	public SortedTreeNode<T> getParent(){
		return parent;
	}
	
	public Set<SortedTreeNode<T>> getChildren(){
		return Collections.unmodifiableSet( children );
	}
	
	@Override
	public int compareTo(SortedTreeNode<T> arg) {
		return data.compareTo(arg.getData());
	}
	
	/*
	 * Returns true if this node is the tree's root
	 */
	public boolean isRoot(){
		return (parent == null);
	}
	
	/*
	 * Returns true if this node is a leaf
	 */
	public boolean isLeaf(){
		return children.isEmpty();
	}
	
	public String toString(){
		String out = "STN (" + data + ") " ;
		
		if( isRoot() ){
			out += "ROOT "; 
		}
		
		// it's possible to be a both a root and a leaf
		// in the case of a tree with a single node
		if( isLeaf() ){
			out += "LEAF ";
		}else{
			out += "[" + children.size() + " children]";
		}
		
		return out;
	}
}