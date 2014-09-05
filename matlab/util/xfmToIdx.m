function idxOut = xfmToIdx( T, sz, rsout )
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

% [ xyz ] = ndgrid ( -half(1) : half(1), ...
%                    -half(2) : half(2), ...
%                    -half(3) : half(3));

docell = 0;
if( iscell(T))
    thisT = T;
    docell = 1;
    idxOut = {};
else
    thisT = { T };
end

npts = numel( xyz{1} );
pts = zeros( ndim, npts );
for i = 1:ndim
   pts( i, : ) = reshape( xyz{i}, 1, [] ); 
end

for n = 1:length( thisT )

    ptXfm = thisT{n} * pts ; 
    ptXfm = ptXfm + repmat( half', 1, size(ptXfm,2)) + 1;

    args = cell( ndim, 1);
    for m = 1:ndim;
        args{m} = ptXfm(m,:);
    end

    if( rsout )
        itmp = reshape( sub2ind( sz, args{:}), sz );
    else
        itmp = sub2ind( sz, args{:} )';
    end


    if( docell )
        idxOut{n} = itmp;
    else
        idxOut = itmp;
    end
end
