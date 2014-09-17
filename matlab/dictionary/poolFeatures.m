function out = poolFeatures( fun, in, num )
% poolFeatures applies a pooling function over a feature matrix
%
% Usage:
%   out = poolFeatureS( in, num );
%
%
% Inputs:
%   fun - a string in {'max', 'mean'} 
%   in  - an array
%         rows represent feature variables 
%         columns represent observations
%         every 'num' columns will be pooled
%   num 
%
% Example -
%   in = rand( 20, 6 );
%   out = poolFeatures( 'max', in, 2 );
%   % out is a 20 x 3 matrix

switch ( fun )
    case 'mean'
        out = reshape( mean( reshape( in', num, [] )', 2 ), size(in,2)./num, [] )';
    case 'max'
        out = reshape( max( reshape( in', num, [] )', [], 2 ), size(in,2)./num, [] )';
    otherwise
        error( 'invalid pooling function option - must be ''mean'' or ''max''');
end

