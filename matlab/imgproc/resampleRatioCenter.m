% function out = resampleRatioCenter( im, ratio )
% % out = resampleRatioCenter( im, ratio )
% % 
% 
% szin  = size( im );
% szout = szin ./ ratio;

%%

rng = -10:10;
[x,y] = meshgrid(rng, rng);

sz = size(x);

dsFactors = [1 2];
radSub = ceil( (sz-1)./ 2 ) ./ dsFactors
del = [ 1 2 ];

[xs, ys] = meshgrid( rng(1):del(1):rng(end), rng(1):del(2):rng(end));