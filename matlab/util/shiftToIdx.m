function [i,sz] = shiftToIdx( shift, subsz, sz )
% [i,sz] = shiftToIdx( shift, subsz, sz )
%
% Generates index lookups for a patch of size 'sz'
% from a shift vector 'shift'
%
% See also:
%   cubeSymmetry
%   xfmToIdx
%
% John Bogovic
% August 2014

% shift should be a column vector

docell = 0;        
if ( iscell( shift ))
    docell = 1;        
    shiftList = shift;
    N = length( shiftList );
elseif( min( size( shift )) > 1)
    shiftList = shift;
    N = size( shiftList, 1 );
else
    shiftList = { shift };
    N = 1;
end

if( ~exist('sz', 'var') || isempty( sz ))
    if( docell )
        maxshift = max( cellfun( @max, shiftList )); 
    else
        maxshift = max( shift );
    end
    sz = subsz + 2*maxshift;
end

if( docell )
    i = {};
%else
%    i = zeros( N, prod(sz) ); 
end

for n = 1:N

    if( docell )
        thisshift = vecToRowCol( shiftList{n}, 'col' );
    else
        thisshift = shiftList(n,:);
    end

    half    = (sz - 1) / 2;
    subhalf = (subsz - 1) / 2;

    num = prod( sz );
    j = reshape( 1:num, sz );
    itmp = circshift( j, thisshift );

    if( N > 1 )
        i{n} = itmp;
    %elseif( N > 1 )
    %    i(n,:) = itmp(:);
    else
        i = itmp;
    end

end
