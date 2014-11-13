function [faceLists, vertLists] = renderPointSpheres( pts, figHandle, ...
                                        rad, res, patchargs )
% Usage:
%   [faceLists, vertLists] = renderPointSpheres( pts, figHandle, rad,
%   res, patchargs )

dorender = 1;
if( ~exist('figHandle','var') || isempty( figHandle ))
    dorender = 1;
end

if( ~exist('res','var') || isempty( res ))
    res = 10;
end
if( ~exist('rad','var') || isempty( rad ))
    rad = 0.5;
end
if( ~exist('patchargs','var') || isempty( patchargs ))
    patchargs = {};
end

% generate vertices and faces for 
[X,Y,Z] = sphere( res );
[sphereF, sphereV] = surf2patch( X, Y, Z );
clear X Y Z;

nV = size( sphereV, 1 );
sv = [sphereV, ones( nV, 1)];

% base transformation helper
Tbase = eye(4);
scales = rad .*[1 1 1];
Tbase(1:3,1:3) = diag( scales );

N = size(pts,1);

faceLists = cell( N, 1 );
vertLists = repmat( {sphereV}, N , 1); 

for i = 1:N

    T = Tbase;
    T( 1:3, 4 ) = pts(i,:);

    v = T * sv';
    v = v(1:3,:);

    if( dorender )
        patch( 'faces', sphereF, 'vertices', v', patchargs{:} );
        hold on;
    end

end

