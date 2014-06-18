function imout = reduceEffectiveResolution( im, factor )
% imout = reduceEffectiveResolution( im, factor )

sz = size( im );
imout = zeros( sz );

for z = 1:factor:sz(3)
   
    zrng = z : (z + factor - 1);
    im_avgz = mean( im(:,:, zrng), 3);
    
    imout(:,:,zrng ) = repmat( im_avgz, [ 1 1 factor ] );
    
end