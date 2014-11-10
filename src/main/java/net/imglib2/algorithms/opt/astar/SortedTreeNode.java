package net.imglib2.algorithms.opt.astar;

import java.io.Serializable;
import java.util.*;


public class SortedTreeNode<T extends Comparable<T>> implements Comparable<SortedTreeNode<T>>,  Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = -6683027733135726639L;
	private T							data;
	private SortedTreeNode<T> 			parent;
	private SortedSet<SortedTreeNode<T>> children;
	
	private int depth; // the depth this node in the tree
					   // root has depth 0 
	
	/**
	 * Constructor
	 * @param data
	 * @param parent
	 */
	public SortedTreeNode(T data, SortedTreeNode<T> parent ){
		this.data = data;
		this.parent = parent;
		
		children = new TreeSet<SortedTreeNode<T>>();
		
		if( parent != null )
			parent.addChild( this );
		
		if( isRoot() )
			depth = 0;
		else
			depth = parent.depth + 1;
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
	
	public boolean removeChild( SortedTreeNode<T> child ){
		return children.remove( child );
	}
	
	public boolean removeChildren( Collection<SortedTreeNode<T>> children ){
		return children.removeAll( children );
	}
	
	public int numDescendants(){
		if( isLeaf() ){
			return 0;
		}else{
			int N = 0;
			Iterator<SortedTreeNode<T>> it = getChildren().iterator();
			while( it.hasNext() ){
				N += it.next().numDescendants();
			}
			return N;
		}
	}
	
	
	public T getData(){
		return data;
	}
	public int getDepth(){
		return depth;
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
	
	public void printAsTree( ){
		String pad = "";
		printAsTree( pad );
	}
	
	public void printAsTree( String pad ){
		
		System.out.println( pad + this );
		
		if( !isLeaf()){
			String nextPad = pad+"\t";
			Iterator<SortedTreeNode<T>> it = getChildren().iterator();
			while( it.hasNext() ){
				it.next().printAsTree(nextPad);
			}
		}
		
	}
	
	public LinkedList<SortedTreeNode<T>> pathToRoot(){
		LinkedList<SortedTreeNode<T>> path = new LinkedList<SortedTreeNode<T>>( );
		
		SortedTreeNode<T> obj = this;
		path.addFirst( obj );
		
		while( !obj.isRoot() ){
			obj = obj.getParent();
			path.addFirst( obj );
		}
		
		return path;
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
	
	public static <T extends Comparable<T>> String pathToString( Collection<SortedTreeNode<T>> path, String pad ){
		String out = "";
		Iterator<SortedTreeNode<T>> it = path.iterator();

		String prefix = "";
		while( it.hasNext() ){
			out = out + prefix + it.next() + "\n";
			prefix = prefix + pad;
		}
		
		return out;
	}
}