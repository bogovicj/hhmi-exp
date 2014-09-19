function [farIdxs, keepIdxs, pd] = dissimilarVecs( X, tol, dist )
% farIdxs = dissimilarVecs( X, tol, dist )
% 
% X    - rows are observations, columns are variables
% tol  - tolerance below which vectors are 
% dist - type of distance function (default = 'euclidean')

if( ~exist( 'tol', 'var') || isempty( tol ))
    tol = 0.1;
end

if( ~exist( 'dist', 'var') || isempty( dist ))
    dist = 'euclidean';
end

pd = pdist( X, dist );
N = size(X,1);
[i,j] = ind2subTri( N, find(pd < tol) );

u = unique( [ i j ] );

farIdxs = [];
for k = u
    if( nnz( farIdxs == k ) > 0 )
       continue; 
    end
    farIdxs = [ farIdxs j(i==k) ];
end

farIdxs = sort( farIdxs );
keepIdxs = setdiff( 1:N, farIdxs );