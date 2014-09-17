function [patches, coords, numNeighbors ] = grabPatchesNeighborhood( im, patchSize, numPatches, coords, mask, nbrhood )
% Usage:
%   [patches, coords] = grabPatchesNeighborhood( im, patchSize, numPatches, coords, mask, nbrhood );

pNumel = prod( patchSize );
pRad = (patchSize - 1)./2;
sz = size( im );

if( ~exist('nbrhood','var') || isempty(nbrhood))
    [x,y,z] = ndgrid(-1:1, -1:1, -1:1);
    nbrhood = [ x(:) y(:) z(:) ];
end

numNeighbors = size( nbrhood, 1 );
% remove 'center' of neighborhood
nbrhood = nbrhood( ...
    setdiff( 1:numNeighbors, ...
             find((nbrhood(:,1)==0) & (nbrhood(:,2)==0) & (nbrhood(:,3)==0))), ...
	:);

numNeighbors = size( nbrhood, 1 );

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

patches = zeros( (numNeighbors + 1).*numPatches, pNumel, class(im));

k = 1;
for i = 1:numPatches
   
    thisPatch = im( x(i)-pRad(1) : x(i)+pRad(1), ...
                    y(i)-pRad(2) : y(i)+pRad(2), ...
                    z(i)-pRad(3) : z(i)+pRad(3) );
                
    patches(k,:) = thisPatch(:);
    k = k + 1;
    
    for j = 1:numNeighbors 
    
        thisPatch = im( x(i)+nbrhood(j,1)-pRad(1) : x(i)+nbrhood(j,1)+pRad(1), ...
                        y(i)+nbrhood(j,2)-pRad(2) : y(i)+nbrhood(j,2)+pRad(2), ...
                        z(i)+nbrhood(j,3)-pRad(3) : z(i)+nbrhood(j,3)+pRad(3) );
                
        patches(k,:) = thisPatch(:);
        k = k + 1;
    end
    
end

coords = {x;y;z};