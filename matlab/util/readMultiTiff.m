function I = readMultiTiff( fn )
% readMultiTiff reads a tiff stack
% Usage I = readMultiTiff( fn )
%
% fn - the file name

i_info = imfinfo(fn);

N = numel(i_info);
tmp = imread(fn);
I = zeros( [ size(tmp) N ]);
for i=1:N
    I(:,:,i) = imread( fn, i, 'Info', i_info );
end