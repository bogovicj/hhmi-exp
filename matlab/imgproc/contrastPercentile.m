function C = contrastPercentile( X, percentileMin, percentileMax, dim )
% CONTRASTPERCENTILE compute image contrast as
%   ( max - min )./ ( max + min )
%   where the max and min are computed using the input percentileMax
%   and percentileMin.  These default to 95 and 5, respectively.
%
% Usage:
% C = contrast( X, percentileMin, percentileMax, dim )
%
% If dim is given the percentiles are computed along the dim^th dimension
% of the input matrix X.  
%   Example:
%       contrast( X, percentileMin, percentileMax, 2 );
%           returns the values of contrast for patches arranged in rows of X 

if( ~exist('dim', 'var'))
    dim = [];
end
if( ~exist('percentileMin', 'var') || isempty( percentileMin ))
    percentileMin = 5;
end
if( ~exist('percentileMax', 'var') || isempty( percentileMax ))
    percentileMax = 95;
end

if( isempty( dim ))
    [ vals ] = prctile( X(:), [percentileMin percentileMax]);
    smallval = vals(1);
    bigval   = vals(2);
else
    [ vals ] = prctile( X, [percentileMin percentileMax], dim );
    if( dim == 1 )
        smallval = vals(1,:);
        bigval   = vals(2,:);
    else
        smallval = vals(:,1);
        bigval   = vals(:,2);
    end
end

C = ( bigval - smallval ) ./ ( bigval + smallval );