package net.imglib2.algorithms.gradient;

import org.ejml.data.DenseMatrix64F;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import ij.gui.Roi;
import ij.io.RoiDecoder;
import edu.jhu.ece.iacl.utility.ArrayUtil;
import net.imglib2.Cursor;
import net.imglib2.ExtendedRandomAccessibleInterval;
import net.imglib2.Localizable;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.io.ImgOpener;
import net.imglib2.type.BooleanType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ImgUtil;
import net.imglib2.util.Intervals;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class PartialDerivativeDirectional {

	public static void main(String[] args){

		long nx = 50;
		long ny = 50;
		
//		ArrayImgFactory<FloatType> factory = new ArrayImgFactory<FloatType>();
//		Img<FloatType> img = genGradient(factory, new FloatType(0), nx, ny, 0.05);
//		
//		Img<FloatType> grad = factory.create(new long[]{nx, ny, img.numDimensions()}, new FloatType(0f));
//		System.out.println("grad: " + grad);

//		testDirDeriv(img, grad);
		
//		testMaskDeriv(img,grad);
//		testMaskDeriv2(img,grad);
		
//		testSingleGradient();
		
		testMaskDerivBin();
		
		System.out.println("finished");
	}
	
	
	public static void testMaskDerivBin(){
		
		long nx = 5;
		long ny = 5;
		
		ArrayImgFactory<BitType> factory = new ArrayImgFactory<BitType>();
		Img<BitType> img = factory.create(new long[]{nx, ny}, new BitType(false));
		
		RandomAccess<BitType> imra = img.randomAccess();
		
		imra.setPosition(new int[]{0,0});
		imra.get().set(true);
		
		imra.setPosition(new int[]{1,0});
		imra.get().set(true);
		
		imra.setPosition(new int[]{0,1});
		imra.get().set(true);
		
		imra.setPosition(new int[]{1,1});
		imra.get().set(true);
		
		imra.setPosition(new int[]{2,2});
		double gval = gradientMask(img, img, 0, imra, false);
		
		System.out.println("gval: " + gval);
		
	}
	
	public static void testSingleGradient(){
		String filename = "/groups/jain/home/bogovicj/learning/advanced-imglib2/images/bee-1.tif";

		UnsignedByteType type = new UnsignedByteType();
		ArrayImgFactory< UnsignedByteType > factory = new ArrayImgFactory< UnsignedByteType >();
		
		Img<UnsignedByteType> img = null;
		Img<BitType>          mask = null;
		Img<FloatType>        gradient = null;
		
		// read image and roi
		try{
			
			img = new ImgOpener().openImg( filename , factory, type );
			mask = img.factory().imgFactory(new BitType(true)).create(img, new BitType(true));
			gradient = img.factory().imgFactory(new FloatType(0f)).create(
					new long[]{img.dimension(0), img.dimension(1), 2}, new FloatType(0f));
		
		}catch(Exception e){
			e.printStackTrace();
		}
		
		System.out.println("grad: " + gradient);
		
		ImgUtil.fill(mask, new BitType(true));
//		ImageJFunctions.show( mask );
		
//		RandomAccessibleInterval< T> infinite1 =
//				Views.interval( Views.extendValue( img,  img.randomAccess().get()), img);
		
		RandomAccessibleInterval<BitType> maskInf = Views.interval(
				Views.extendValue(mask, new BitType(false)), img);
		
		RandomAccess<BitType> mcra = mask.randomAccess();
		mcra.setPosition(new int[]{25,25});
		
		
		
//		while(mc.hasNext()){
		
			gradientMask(img, maskInf, gradient, 0, mcra, true);
			gradientMask(img, maskInf, gradient, 1, mcra, true);
//		}
			
		
		
		ImgUtil.printNumNonZero(gradient);
		
//		System.out.println(" ");
//		
//		RandomAccess<FloatType> gdra = gradient.randomAccess();
//		
//		gdra.setPosition(new int[]{25,25,0});
//		System.out.println(" " + gdra.get());
//		
//		gdra.setPosition(new int[]{25,25,1});
//		System.out.println(" " + gdra.get());
	}
	
	public static Img<FloatType> genGradient( ArrayImgFactory<FloatType> factory, FloatType val,  long nx, long ny, double stdDev){
		
		Img<FloatType> img = factory.create(new long[]{nx, ny}, val);
		
		boolean a = true;
		RandomAccess<FloatType> ra = img.randomAccess();
		for(long x=0; x<nx; x++){
			for(long y=0; y<ny; y++){
		
			ra.setPosition(new long[]{x, y});
			if(a){
				ra.get().set((float)x + 1.2f + (float)(stdDev*Math.random()));
			}else{
				ra.get().set((float)x + 0.8f + (float)(stdDev*Math.random()));	
			}
			
			}
			a = !a;
		}
		
		return img;
	}
	
	public static <T extends RealType<T>> void testMaskDeriv2 ( Img<T> img, Img<T> grad ){
		ArrayImgFactory<BitType> factory = new ArrayImgFactory<BitType>();
		long nx = img.dimension(0);
		long ny = img.dimension(1);
		
		int fixedy = 1;
		int cx     = 1;
		Img<BitType> mask = factory.create(new long[]{nx, ny}, new BitType(true));
		ImgUtil.fill(mask, new BitType(true));
		RandomAccess<BitType> maskRa = mask.randomAccess();
		maskRa.setPosition(new int[]{cx, fixedy});
		
		int d = 0;
		PartialDerivativeDirectional.gradientMask(img, mask, grad, d, maskRa, true);
	}
	
	public static <T extends RealType<T>> void testMaskDeriv ( Img<T> img, Img<T> grad ){
		ArrayImgFactory<BitType> factory = new ArrayImgFactory<BitType>();
		long nx = img.dimension(0);
		long ny = img.dimension(1);
		
		int fixedy = 1;
		int cx     = 1;
		Img<BitType> mask = factory.create(new long[]{nx, ny}, new BitType(true));
		ImgUtil.fill(mask, new BitType(true));
		RandomAccess<BitType> maskRa = mask.randomAccess();
		
		Cursor<BitType> maskCursor = mask.cursor();
		while(maskCursor.hasNext()){
			if(maskCursor.next().get()){
				System.out.println("true!");
				break;
			}
		}
		
		System.out.println(mask.numDimensions());
//		System.out.println(mask.dimension(0));
//		System.out.println(mask.dimension(1));
		
//		maskRa.setPosition(new int[]{cx-1,fixedy});
//		maskRa.get().set(new BitType(false));
		
//		maskRa.setPosition(new int[]{cx+1,fixedy});
//		maskRa.get().set(new BitType(false));
		
		RandomAccess<T> raGrad = grad.randomAccess();
		raGrad.setPosition(new int[]{cx,fixedy, 0});
		
		maskRa.setPosition(new int[]{cx,fixedy});
		
		gradientMask(img, mask, grad, 0, maskRa, true);
		
		System.out.println("val: " + raGrad.get());
		
	}

	public static <T extends NumericType<T>> void testDirDeriv(Img<T> img, Img<T> grad){
		//	    ExtendedRandomAccessibleInterval view = Views.extendMirrorSingle(img);

		//		gradientCentralDifference(view, Views.hyperSlice(grad, 2, 0), 0);
		//		gradientBackwardDifference(view, Views.hyperSlice(grad, 2, 0), 0);
		//		gradientForwardDifference(view, Views.hyperSlice(grad, 2, 0), 0);

		//		gradientCentralDifference(view, Views.hyperSlice(grad, 2, 1), 1);


		RandomAccess<T> raImg = img.randomAccess();
		raImg.setPosition(new int[]{25,25});

		RandomAccess<T> raGrad = grad.randomAccess();
		raGrad.setPosition(new int[]{25,25,0});


		//		gradientForwardDifference(img, grad, 0, raGrad);
		//		gradientForwardDifference(img, grad, 1, raGrad);

		gradientBackwardDifference(img, grad, 0, raGrad);
		gradientBackwardDifference(img, grad, 1, raGrad);

		//		raGrad = grad.randomAccess();
		//		raGrad.setPosition(new int[]{25,25,0});
		System.out.println("\n " + raGrad.get());
		//		raGrad.fwd(2);
		//		System.out.println(" " + raGrad.get());



		//		ImageJFunctions.show( img );
		//		ImageJFunctions.show(  Views.hyperSlice(grad, 2, 0) );
		//		ImageJFunctions.show(  Views.hyperSlice(grad, 2, 1) );

		//		IntervalView<FloatType> myview = Views.hyperSlice( Views.hyperSlice(grad, 2, 0), 0, 25);
		////		IntervalView<FloatType> myview = Views.hyperSlice( img, 0, 0);
		//		
		//		RandomAccess<FloatType> rasamp = myview.getSource().randomAccess();
		//		System.out.println("rs ndims: " + rasamp.numDimensions());
		//		
		//		for(int i=0; i<20; i++){
		//			rasamp.fwd(0);
		//			System.out.println(" " + rasamp.get().get());
		//		}
	}
	

	// fast version
	/**
	 * Compute the partial derivative of source in a particular dimension.
	 *
	 * @param source
	 *            source image, has to provide valid data in the interval of the
	 *            gradient image plus a one pixel border in dimension.
	 * @param gradient
	 *            output image
	 * @param dimension
	 *            along which dimension the partial derivatives are computed
	 */
	public static < T extends NumericType< T > > void gradientCentralDifference( final RandomAccessible< T > source, final RandomAccessibleInterval< T > gradient, final int dimension )
	{
		final int n = gradient.numDimensions();

		final long[] min = new long[ n ];
		gradient.min( min );
		final long[] max = new long[ n ];
		gradient.max( max );
		final long[] shiftback = new long[ n ];
		for ( int d = 0; d < n; ++d )
			shiftback[ d ] = min[ d ] - max[ d ];

		final RandomAccess< T > result = gradient.randomAccess();
		final RandomAccess< T > back = source.randomAccess( Intervals.translate( gradient, 1, dimension ) );
		final RandomAccess< T > front = source.randomAccess( Intervals.translate( gradient, -1, dimension ) );

		result.setPosition( min );
		back.setPosition( min );
		back.bck( dimension );
		front.setPosition( min );
		front.fwd( dimension );

		final long max0 = max[ 0 ];
		while ( true )
		{
			// process pixel
			final T t = result.get();
			t.set( front.get() );
			t.sub( back.get() );
			t.mul( 0.5 );

			// move to next pixel
			// check dimension 0 separately to avoid the loop over d in most iterations
			if ( result.getLongPosition( 0 ) == max0 )
			{
				if ( n == 1 )
					return;
				result.move( shiftback[ 0 ], 0 );
				back.move( shiftback[ 0 ], 0 );
				front.move( shiftback[ 0 ], 0 );
				// now check the remaining dimensions
				for ( int d = 1; d < n; ++d )
					if ( result.getLongPosition( d ) == max[ d ] )
					{
						result.move( shiftback[ d ], d );
						back.move( shiftback[ d ], d );
						front.move( shiftback[ d ], d );
						if ( d == n - 1 )
							return;
					}
					else
					{
						result.fwd( d );
						back.fwd( d );
						front.fwd( d );
						break;
					}
			}
			else
			{
				result.fwd( 0 );
				back.fwd( 0 );
				front.fwd( 0 );
			}
		}
	}

	// fast version
	/**
	 * Compute the partial derivative of source in a particular dimension.
	 *
	 * @param source
	 *            source image, has to provide valid data in the interval of the
	 *            gradient image plus a one pixel border in dimension.
	 * @param gradient
	 *            output image
	 * @param dimension
	 *            along which dimension the partial derivatives are computed
	 */
	public static < T extends NumericType< T > > void gradientBackwardDifference( final RandomAccessible< T > source, final RandomAccessibleInterval< T > gradient, final int dimension )
	{
		final int n = gradient.numDimensions();

		final long[] min = new long[ n ];
		gradient.min( min );
		final long[] max = new long[ n ];
		gradient.max( max );
		final long[] shiftback = new long[ n ];
		for ( int d = 0; d < n; ++d )
			shiftback[ d ] = min[ d ] - max[ d ];

		final RandomAccess< T > result = gradient.randomAccess();
		final RandomAccess< T > back = source.randomAccess( Intervals.translate( gradient, 1, dimension ) );
		final RandomAccess< T > ctr = source.randomAccess( gradient );

		result.setPosition( min );
		back.setPosition( min );
		back.bck( dimension );
		ctr.setPosition( min );

		final long max0 = max[ 0 ];
		while ( true )
		{
			// process pixel
			final T t = result.get();
			t.set( ctr.get() );
			t.sub( back.get() );

			// move to next pixel
			// check dimension 0 separately to avoid the loop over d in most iterations
			if ( result.getLongPosition( 0 ) == max0 )
			{
				if ( n == 1 )
					return;
				result.move( shiftback[ 0 ], 0 );
				back.move( shiftback[ 0 ], 0 );
				ctr.move( shiftback[ 0 ], 0 );
				// now check the remaining dimensions
				for ( int d = 1; d < n; ++d )
					if ( result.getLongPosition( d ) == max[ d ] )
					{
						result.move( shiftback[ d ], d );
						back.move( shiftback[ d ], d );
						ctr.move( shiftback[ d ], d );
						if ( d == n - 1 )
							return;
					}
					else
					{
						result.fwd( d );
						back.fwd( d );
						ctr.fwd( d );
						break;
					}
			}
			else
			{
				result.fwd( 0 );
				back.fwd( 0 );
				ctr.fwd( 0 );
			}
		}
	}
	
	public static < T extends NumericType< T > > void gradientBackwardDifference( final RandomAccessible< T > source, final RandomAccessibleInterval< T > gradient, final int dimension, Localizable pos )
	{
		final int n = gradient.numDimensions();

		final long[] min = new long[ n ];
		for ( int d = 0; d < n; ++d ){
			if(d < (n-1)){
				min[d] = pos.getLongPosition(d);
			}else{
				min[d] = dimension;
			}
		}

		final RandomAccess< T > result = gradient.randomAccess();
		final RandomAccess< T > back = source.randomAccess( Intervals.translate( gradient, 1, dimension ) );
		final RandomAccess< T > ctr = source.randomAccess( gradient );

		result.setPosition( min );
		back.setPosition( min );
		back.bck( dimension );
		ctr.setPosition( min );

		// process pixel
		final T t = result.get();
		t.set( ctr.get() );
		t.sub( back.get() );

//		System.out.println(" bd " + result.get());
	}
	
	// fast version
	/**
	 * Compute the partial derivative of source in a particular dimension using a
	 * forward difference scheme.
	 *
	 * @param source
	 *            source image, has to provide valid data in the interval of the
	 *            gradient image plus a one pixel border in dimension.
	 * @param gradient
	 *            output image
	 * @param dimension
	 *            along which dimension the partial derivatives are computed
	 */
	public static < T extends NumericType< T > > void gradientForwardDifference( final RandomAccessible< T > source, final RandomAccessibleInterval< T > gradient, final int dimension )
	{
		final int n = gradient.numDimensions();

		final long[] min = new long[ n ];
		gradient.min( min );
		final long[] max = new long[ n ];
		gradient.max( max );
		final long[] shiftback = new long[ n ];
		for ( int d = 0; d < n; ++d )
			shiftback[ d ] = min[ d ] - max[ d ];

		final RandomAccess< T > result = gradient.randomAccess();
		final RandomAccess< T > ctr = source.randomAccess(  gradient );
		final RandomAccess< T > front = source.randomAccess( Intervals.translate( gradient, -1, dimension ) );

		result.setPosition( min );
		
		ctr.setPosition( min );
		
		front.setPosition( min );
		front.fwd( dimension );

		final long max0 = max[ 0 ];
		while ( true )
		{
			// process pixel
			final T t = result.get();
			t.set( front.get() );
			t.sub( ctr.get() );

			// move to next pixel
			// check dimension 0 separately to avoid the loop over d in most iterations
			if ( result.getLongPosition( 0 ) == max0 )
			{
				if ( n == 1 )
					return;
				result.move( shiftback[ 0 ], 0 );
				ctr.move( shiftback[ 0 ], 0 );
				front.move( shiftback[ 0 ], 0 );
				// now check the remaining dimensions
				for ( int d = 1; d < n; ++d )
					if ( result.getLongPosition( d ) == max[ d ] )
					{
						result.move( shiftback[ d ], d );
						ctr.move( shiftback[ d ], d );
						front.move( shiftback[ d ], d );
						if ( d == n - 1 )
							return;
					}
					else
					{
						result.fwd( d );
						ctr.fwd( d );
						front.fwd( d );
						break;
					}
			}
			else
			{
				result.fwd( 0 );
				ctr.fwd( 0 );
				front.fwd( 0 );
			}
		}
	}

	// fast version
		/**
		 * Compute the partial derivative of source in a particular dimension using a
		 * forward difference scheme.
		 *
		 * @param source
		 *            source image, has to provide valid data in the interval of the
		 *            gradient image plus a one pixel border in dimension.
		 * @param gradient
		 *            output image
		 * @param dimension
		 *            along which dimension the partial derivatives are computed
		 */
		public static < T extends NumericType< T > > void gradientForwardDifference( final RandomAccessible< T > source, final RandomAccessibleInterval< T > gradient, final int dimension, Localizable pos )
		{
			final int n = gradient.numDimensions();

			final long[] min = new long[ n ];
//			System.out.println("pos: " );
			for ( int d = 0; d < n; ++d ){
				if(d < (n-1)){
					min[d] = pos.getLongPosition(d);
				}else{
					min[d] = dimension;
				}
//				System.out.println("  " + min[d]);
			}
//			System.out.println("  ");


			final RandomAccess< T > result = gradient.randomAccess();
			final RandomAccess< T > ctr = source.randomAccess(  gradient );
			final RandomAccess< T > front = source.randomAccess( Intervals.translate( gradient, -1, dimension ) );

			result.setPosition( min );
			
			ctr.setPosition( min );
			
			front.setPosition( min );
			front.fwd( dimension );

			// process pixel
			final T t = result.get();
			t.set( front.get() );
			t.sub( ctr.get() );

//			System.out.println(" fd " + result.get());
			
		}
		
		/**
		 * Compute the partial derivative of source in a particular dimension using a
		 * forward difference scheme.
		 *
		 * @param source
		 *            source image, has to provide valid data in the interval of the
		 *            gradient image plus a one pixel border in dimension.
		 * @param gradient
		 *            output image
		 * @param dimension
		 *            along which dimension the partial derivatives are computed
		 */
		public static < T extends RealType< T >, S extends RealType< S > , B extends BooleanType<B>> void gradientMask( final RandomAccessible< T > source, final RandomAccessible<B> mask, final RandomAccessibleInterval< S > gradient, final int dimension, Localizable pos, boolean maskIsTrue )
		{
			int n = gradient.numDimensions();
			
			final RandomAccess< S > result = gradient.randomAccess();
			final RandomAccess< T > srcRA = source.randomAccess( gradient );
			final RandomAccess< B > mskRA = mask.  randomAccess( gradient );

			result.setPosition( pos ); // res at 0
			result.setPosition( dimension, (n-1) );
			
			srcRA.setPosition( pos );  // src at 0
			mskRA.setPosition( pos );  // msk at 0
			
//			printPosition(result);
//			printPosition(srcRA);
//			printPosition(mskRA);
			boolean maskVal = optNegateBoolean( mskRA.get().get(), maskIsTrue);
			
			
			if( ! maskVal ){
				return;
			}
			
			// there's something to do
			final S r = result.get();
			r.setZero();
			
			// check fwd deriv
			boolean didFwd = false;
			mskRA.fwd(dimension); // msk at 1
			
			maskVal = optNegateBoolean( mskRA.get().get(), maskIsTrue);
			
			if( maskVal ){
//				System.out.println("doing fwd");
				
				double d0 = srcRA.get().getRealDouble();
				srcRA.fwd(dimension);  // src at 1
				double d1 = srcRA.get().getRealDouble();
				
				
				r.setReal(d1-d0);
				
				didFwd = true;
			}
			
			srcRA.setPosition( pos ); // src at 0
			
			mskRA.bck(dimension); // msk at 0
			mskRA.bck(dimension); // msk at -1
			
			maskVal = optNegateBoolean( mskRA.get().get(), maskIsTrue);
			
			if( maskVal ){
//				System.out.println("doing rev");
				
				double d0 = srcRA.get().getRealDouble();
				srcRA.bck(dimension); // src at -1
				double dm1 = srcRA.get().getRealDouble();
				
				
				if(didFwd){ // need to normalize
					r.setReal( 0.5 * (r.getRealDouble() + d0 - dm1) );
				}else{
					r.setReal(d0 - dm1);
				}
			}
			
//			System.out.println("grad val real double: " + r.getRealDouble());
//			System.out.println("grad val type: " + r.copy());
			
		}
		
		/**
		 * Compute the partial derivative of source in a particular dimension using a
		 * forward difference scheme.
		 *
		 * @param source
		 *            source image, has to provide valid data in the interval of the
		 *            gradient image plus a one pixel border in dimension.
		 * @param gradient
		 *            output image
		 * @param dimension
		 *            along which dimension the partial derivatives are computed
		 */
		public static < T extends RealType< T >, S extends RealType< S > , B extends BooleanType<B>> double gradientMask( final RandomAccessible< T > source, final RandomAccessible<B> mask, final int dimension, Localizable pos, boolean maskIsTrue )
		{
			
			final RandomAccess< T > srcRA = source.randomAccess( );
			final RandomAccess< B > mskRA = mask.  randomAccess( );

			
			srcRA.setPosition( pos );  // src at 0
			mskRA.setPosition( pos );  // msk at 0
			
//			printPosition(result);
//			printPosition(srcRA);
//			printPosition(mskRA);
			
			boolean maskVal = optNegateBoolean( mskRA.get().get(), maskIsTrue);
			if( ! maskVal ){
				return Double.NaN;
			}
			
			// there's something to do
			double r = 0;
			
			// check fwd deriv
			boolean didFwd = false;
			
			mskRA.fwd(dimension); // msk at 1
			maskVal = optNegateBoolean( mskRA.get().get(), maskIsTrue);
			
			if( maskVal ){
//				System.out.println("doing fwd");
				
				double d0 = srcRA.get().getRealDouble();
				srcRA.fwd(dimension);  // src at 1
				double d1 = srcRA.get().getRealDouble();
				
				r = d1 - d0;
				
				didFwd = true;
			}
			
			srcRA.setPosition( pos ); // src at 0
			
			mskRA.bck(dimension); // msk at 0
			mskRA.bck(dimension); // msk at -1
			
			maskVal = optNegateBoolean( mskRA.get().get(), maskIsTrue);
			
			if( maskVal ){
//				System.out.println("doing rev");
				
				double d0 = srcRA.get().getRealDouble();
				srcRA.bck(dimension); // src at -1
				double dm1 = srcRA.get().getRealDouble();
				
				if(didFwd){ // need to normalize
					r = ( 0.5 * (r + d0 - dm1) );
				}else{
					r = d0 - dm1;
				}
			}
			
//			System.out.println("r: " + r);
			
			return r;
		}
		
		
		public static boolean optNegateBoolean(boolean in, boolean thruSwitch){
			if( thruSwitch ){
				return in;
			}else{
				return !in;
			}
		}
		
		public static void printPosition(RandomAccess ra){
			int N = ra.numDimensions();
			for(int d=0; d<N; d++){
				System.out.print( " " + ra.getLongPosition(d));
			}
			System.out.print("\n");
		}
	
}
