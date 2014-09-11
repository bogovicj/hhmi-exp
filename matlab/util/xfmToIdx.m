function i = xfmToIdx( T, sz, rsout )
% i = xfmToIdx( T, sz, rsout )
%
% Generates index lookups for a patch of size 'sz'
% from a transformation matrix 'T'. 
%
% rsout - "reshape the output" (Default = true)
%
% See also:
%   cubeSymmetry
%
% John Bogovic
% August 2014

if( ~exist( 'rsout', 'var') || isempty( rsout ))
    rsout = 1;
end

ndim = length( sz );
half = (sz - 1) / 2;
xyz = cell( ndim, 1 );
args = cell( ndim, 1 );
for idxOut = 1:ndim
    args{idxOut} = -half(idxOut) : half(idxOut);
end

[ xyz{:} ] = ndgrid ( args{:} );


docell = 0;
if( iscell(T))
    thisT = T;
    docell = 1;
    i = {};
else
    thisT = { T };
end

npts = numel( xyz{1} );
pts = zeros( ndim, npts );
for j = 1:ndim
   pts( j, : ) = reshape( xyz{j}, 1, [] ); 
end

for n = 1:length( thisT )

    % the output size may be different after
    % accountoing for permutation of the dimensions
    szxfm = abs(thisT{n} * sz');
    halfxfm = (szxfm - 1) / 2;
%     perm = abs( thisT{n} * [1:length(sz)]' );
    
    
    ptXfm = thisT{n} * pts; 
    ptXfm = ptXfm + repmat( halfxfm, 1, size(ptXfm,2)) + 1;

    
    args = {};
    for m = 1:size(ptXfm,1)
        args{m} = ptXfm(m,:);
    end
    %args = args( perm) ;

    if( rsout )
        itmp = reshape( sub2ind( szxfm', args{:}), szxfm' );
    else
        itmp = sub2ind( szxfm', args{:} )';
    end


    if( docell )
        i{n} = itmp;
    else
        i = itmp;
    end
end
