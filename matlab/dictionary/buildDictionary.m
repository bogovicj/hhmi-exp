function D = buildDictionary( X, method, params )
% buildDictionary from data
% Usage:
%   D = buildDictionary( X, method, params )
%
%   X - ( N x M )
%       N - number of observations
%       M - number of variables
%
% Rows of the ouput matrix contain 
% 
% methods:
%   spams
%   kmeans
%   pca     : Note, assumes data are already centered
%
% See also KMEANS, PCA, MEXTRAINDL

if( ~exist( 'params', 'var' ))
    params = {};
end

switch( method )
    
    case 'spams'
        D = mexTrainDL( X', params );
        D = D';
    case 'kmeans'
        [~,D] = kmeans( X, params );
    case 'pca'
        D = pca( X, 'centered', false );
        D = D';
    otherwise
        error('Invalid method')
end