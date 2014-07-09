function [ii] = visDictionary( D, patchSize, dir, num )
% visDictionary
%
% ii = visDictionary( D, dir, num )


% dir = '/nobackup/saalfeld/john/dictionary_vis/iso_1k_15';
if( ~exist( dir, 'dir'))
    mkdir(dir);
end
ii = randi(size(D,2), 1, num);

for n = 1:num
    n
    im = reshape( D(:,ii(n)), patchSize);
    imdisp( permute( im, [1 2 4 3]), [min(im(:)), max(im(:))], ...
            'Size', 5, 'Border', 0.01);
    
    fn = sprintf('%s/dictPatch_%04d', dir, ii(n));
    
    export_fig( fn, '-m2' );
    close all;
    
end