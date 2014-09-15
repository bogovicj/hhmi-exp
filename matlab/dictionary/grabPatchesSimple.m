function [patches, coords] = grabPatchesSimple( im, patchSize, numPatches, coords, mask )
% Usage:
%   [patches, coords] = grabPatchesSimple( im, patchSize, numPatches, coords, mask );

pNumel = prod( patchSize );
pRad = (patchSize - 1)./2;
sz = size( im );

if( ~exist('coords','var') || isempty(coords))
    if( ~exist('mask','var') || isempty( mask ))
        samplerng = [ (pRad+1)' (sz - pRad - 1)'];
        x = randi( samplerng(1,:), numPatches, 1);
        y = randi( samplerng(2,:), numPatches, 1);
        z = randi( samplerng(3,:), numPatches, 1);
    else
        j = find( mask );
        [x,y,z] = ind2sub( size(mask), j(randi(length(j), numPatches, 1)));
    end
else
    x = coords{1};
    y = coords{2};
    z = coords{3};
    numPatches = length(x);
end

patches = zeros( numPatches, pNumel, class(im));

for i = 1:numPatches
   
    thisPatch = im( x(i)-pRad(1) : x(i)+pRad(1), ...
                    y(i)-pRad(2) : y(i)+pRad(2), ...
                    z(i)-pRad(3) : z(i)+pRad(3) );
                
    patches(i,:) = thisPatch(:);
    
end

coords = {x;y;z};