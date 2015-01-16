function im = imageFromPatches( pMtx, imsz, patchSize, coords, centered)
% Usage:
%   im = imageFromPatches( pMtx, imsz, patchSize, coords, centered) 

if( ~exist('centered','var') || isempty( centered))
    centered = 1;
end

N = size( pMtx, 1);
im = zeros( imsz );
pRad = (patchSize - 1)./2;

xp = coords{1};
yp = coords{2};
zp = coords{3};
    
for i = 1:N

    if( centered )
        x = xp(i)-pRad(1) : xp(i)+pRad(1);
        y = yp(i)-pRad(2) : yp(i)+pRad(2);
        z = zp(i)-pRad(3) : zp(i)+pRad(3); 
    else
        x = xp(i) : xp(i) + patchSize(1) - 1;
        y = yp(i) : yp(i) + patchSize(2) - 1;
        z = zp(i) : zp(i) + patchSize(3) - 1;
    end

    im( x, y, z ) =  reshape(pMtx(i,:),patchSize);

end
