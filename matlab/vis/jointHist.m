function [histCount,C] = jointHist( dat, N, cmap, make_it_tight )
% jointHist( histCount )
%
% 

if( ~exist('make_it_tight','var') || isempty( make_it_tight ))
    make_it_tight = true;    
end

if( ~exist('cmap','var') || isempty( cmap ))
    cmap = diverging_map(linspace(0,1,256), [0.2 0.2 0.7], [0.7 0.2 0.2]);
end

subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end


maxval = max(dat(:));
minval = min(dat(:));

rng = linspace( minval, maxval, N );

[histCount,C] = hist3( dat, { rng, rng } );

figure('color', 'w' );

subplot(2,2,2);
area( rng, sum(histCount,1) ); axis square;
axtmp = axis;
axis( [minval maxval axtmp(3) axtmp(4)]);
xlabel('var 1');
ylabel('Count');

subplot(2,2,3);
area( sum(histCount,2), rng); axis square;
set( gca, 'xdir', 'reverse');
axtmp = axis;
axis( [axtmp(1) axtmp(2) minval maxval  ]);
xlabel('Count');

subplot(2,2,4);
sc(histCount, cmap); axis square;
axis xy;
ylabel('var 2');