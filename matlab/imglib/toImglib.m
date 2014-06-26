function imglibimg = toImglib( img )
% Returns an net.imglib2.img.array.ArrayImg 
%   imglibimg = toImageLib( img )


sz = int64( size( img ));

imglibimg = net.imglib2.img.array.ArrayImgs.floats( ...
                reshape( img, 1, []), sz );