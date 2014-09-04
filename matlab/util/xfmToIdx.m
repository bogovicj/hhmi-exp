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

half = (sz - 1) / 2;

[x,y,z] = ndgrid ( -half(1) : half(1), ...
                   -half(2) : half(2), ...
                   -half(3) : half(3));

docell = 0;
if( iscell(T))
    thisT = T;
    docell = 1;
    i = {};
else
    thisT = { T };
end

for n = 1:length( thisT )

    ptXfm = thisT{n} * [ x(:)'; y(:)'; z(:)' ]; 
    ptXfm = ptXfm + repmat( half', 1, size(ptXfm,2)) + 1;

    args = {};
    for m = 1:size(ptXfm,1)
        args{m} = ptXfm(m,:);
    end

    if( rsout )
        itmp = reshape( sub2ind( sz, args{:}), sz );
    else
        itmp = sub2ind( sz, args{:} )';
    end


    if( docell )
        i{n} = itmp;
    else
        i = itmp;
    end
end
