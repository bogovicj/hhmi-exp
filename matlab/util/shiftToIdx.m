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
else
    shiftList = { shift };
end

if( ~exist('sz', 'var') || isempty( sz ))
    if( docell )
        maxshift = max( cellfun( @max, shiftList )); 
    else
        maxshift = max( shift );
    end
    sz = subsz + 2*maxshift;
end

for n = 1:length(shiftList) 

    thisshift = vecToRowCol( shiftList{n}, 'col' );


    half    = (sz - 1) / 2;
    subhalf = (subsz - 1) / 2;

    num = prod( sz );
    j = reshape( 1:num, sz );
    itmp = circshift( j, thisshift );

    if( docell )
        i{n} = itmp
    else
        i = itmp;
    end

end
