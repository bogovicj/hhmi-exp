function [patches, coords] = grabPatchesSimple( im, patchSize, numPatches, coords, mask )
% Usage:
%   [patches, coords] = grabPatchesSimple( im, patchSize, numPatches, coords, mask );

if( length( patchSize) == 2 )
    patchSize = [ patchSize 1 ];
end
pNumel = prod( patchSize );
pRad = (patchSize - 1)./2;

sz = size( im );

if( ~exist('coords','var') || isempty(coords))
    
    if( numPatches < 0 )
        fprintf('grabPatchesSimple - choosing all coordinates\n');
        if( ~exist('mask','var') || isempty( mask ))
            mask = true( size( im ));
            if( pRad(1) > 0 )
                mask( 1:pRad(1)+1, :, :) = 0;
                mask( end-pRad(1):end, :, :) = 0;
            end
            if( pRad(2) > 0 )
                mask( :, 1:pRad(2)+1, :) = 0;
                mask( :, end-pRad(2):end, :) = 0;
            end
            if( pRad(3) > 0 )
                mask( :, :, 1:pRad(3)+1) = 0;
                mask( :, :, end-pRad(3):end) = 0;
            end
        end
        [x,y,z] = ind2sub( size(mask), find( mask ));
        numPatches = length(x);
    else
        fprintf('grabPatchesSimple - randomly selecting coordinates\n');
        if( ~exist('mask','var') || isempty( mask ))
            samplerng = [ (pRad+1)' (sz - pRad - 1)'];
            x = randi( samplerng(1,:), numPatches, 1);
            y = randi( samplerng(2,:), numPatches, 1);
            z = randi( samplerng(3,:), numPatches, 1);
        else
            if( pRad(1) > 1 )
                mask( 1:pRad(1)+1, :, :) = 0;
                mask( end-pRad(1):end, :, :) = 0;
            end
            if( pRad(2) > 1 )
                mask( :, 1:pRad(2)+1, :) = 0;
                mask( :, end-pRad(2):end, :) = 0;
            end
            if( pRad(3) > 1 )
                mask( :, :, 1:pRad(3)+1) = 0;
                mask( :, :, end-pRad(3):end) = 0;
            end
            
            j = find( mask );
            [x,y,z] = ind2sub( size(mask), j(randi(length(j), numPatches, 1)));
        end
    end
else
    fprintf('grabPatchesSimple - using input coordinates\n');
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