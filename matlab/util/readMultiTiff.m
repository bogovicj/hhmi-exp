function I = readMultiTiff( fn, zrng )
% readMultiTiff reads a tiff stack
% Usage I = readMultiTiff( fn, zrng )
%
% fn   - the file name
% zrng - the zrange to read (optional, defaults to all available slices)

if( ~exist('zrng','var'))
    zrng = [];
end

i_info = imfinfo(fn);

if( isempty( zrng ))
    N = numel(i_info);
else
    N = length( zrng );
end

tmp = imread(fn);

I = zeros( [ size(tmp) N ]);

if( length(i_info(1).BitsPerSample) == 3) % RGB
    fprintf('RGB\n');
    I = permute( I, [1 2 4 3]);
end

for i=1:N
    if( isempty( zrng ))
        I(:,:,i,:) = imread( fn, i, 'Info', i_info );
    else
        I(:,:,i,:) = imread( fn, zrng(i), 'Info', i_info );
    end
    
end