package net.imglib2.algorithms.boundary;

import net.imglib2.AbstractCursorInt;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.algorithm.region.localneighborhood.Neighborhood;
import net.imglib2.algorithm.region.localneighborhood.Shape;
import net.imglib2.img.Img;
import net.imglib2.type.BooleanType;

/**
 * 
 * This cursor iterates over image locations that are on the boundary of a
 * BooleanType image.  Pixels on the boundary are defined as those with
 * a false value but having a neighboring pixel with a true value.
 * 
 * The neighborhood is defined by the input {@link Shape}.
 * 
 * @author John Bogovic 
 *
 * @param <B>
 */
public class BoundaryCursor< B extends BooleanType<B>> extends AbstractCursorInt< B >{

	Shape 		shape;
	Img<B> 		mask;
	Cursor<B> 	maskCursor;
	
	RandomAccess<Neighborhood<B>> nbrhoodAccess;
	
	B current;
	int jumpsToNextPos;
	boolean hasNext = false;
	
	public BoundaryCursor(int n) {
		super(n);
	}

	public BoundaryCursor(Shape shape, Img<B> mask){
		super(mask.numDimensions());
		this.shape = shape;
		this.mask = mask;
		this.maskCursor = mask.cursor();
		
		nbrhoodAccess = shape.neighborhoodsRandomAccessible(mask).randomAccess();
	}

	
	public BoundaryCursor(int n, Shape shape, Img<B> mask){
		super(n);
		this.shape = shape;
		this.mask = mask;
		this.maskCursor = mask.cursor();
		
		nbrhoodAccess = shape.neighborhoodsRandomAccessible(mask).randomAccess();
	}

	public B get() {
		return maskCursor.get();
	}

	public void fwd() {
		maskCursor.jumpFwd(jumpsToNextPos);
	}

	public boolean hasNext() {
		
		int jumps_tmp = 0;
		int[] posin = new int[mask.numDimensions()];
		
		Cursor<B> tmpCursor = maskCursor.copyCursor();
		
		while(tmpCursor.hasNext() ){
			
			tmpCursor.fwd();
			jumps_tmp++;
			
			
			if(! tmpCursor.get().get()){
				
				nbrhoodAccess.setPosition(tmpCursor);
				
				Neighborhood<B> nbr = nbrhoodAccess.get();		
				
				// iterate over neighborhood
				Cursor<B> nbrCursor = nbr.cursor();
				nbr: while (nbrCursor.hasNext()){ // loop over neighbors

					nbrCursor.fwd();
					nbrCursor.localize(posin);
					
					for(int d=0; d<posin.length; d++){ // skip if neighbor is out-of-bounds 
						if( (posin[d] < 0) || (posin[d] >= mask.dimension(d)) ){ 
							continue nbr; 
						}
					}
					
					if( nbrCursor.get().get() ){	// at a boundary pixel
						
						jumpsToNextPos = jumps_tmp;
						return true;
						
					}
				}
			}
			
		}
		return false;
	}
	
//	@Override
	public boolean hasNextOld() {
		int jumps_tmp = 0;
		while(maskCursor.hasNext()){

			maskCursor.fwd();
			jumps_tmp++;
			
			if(maskCursor.get().get()){ // in mask

				nbrhoodAccess.setPosition(maskCursor);
				Neighborhood<B> nbr = nbrhoodAccess.get();

				// iterate over neighborhood
				Cursor<B> nbrCursor = nbr.cursor();
				while (nbrCursor.hasNext()){

					nbrCursor.fwd();

					if( ! nbrCursor.get().get() ){	// at a boundary pixel
						jumpsToNextPos = jumps_tmp;
						return true;
					}
				}
			}
		}
		return false;
	}

	public void reset() {
		maskCursor.reset();
	}

	public int getIntPosition(int arg0) {
		return maskCursor.getIntPosition(arg0);
	}

	public void localize(int[] arg0) {
		maskCursor.localize(arg0);
	}

	@Override
	public AbstractCursorInt<B> copy() {
		return new BoundaryCursor<B>(super.n, shape, mask);
	}

	@Override
	public AbstractCursorInt<B> copyCursor() {
		return copy();
	}
	

}
