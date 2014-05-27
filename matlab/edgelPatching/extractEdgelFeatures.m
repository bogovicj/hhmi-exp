% exploring matching features for crack edgel patches
%
% run_script('extractEdgelFeatures', 'FREAK features ds 4');
% run_script('extractEdgelFeatures', 'FREAK features ds 4, testPt');

%% setup
global SAVEPATH;

ds = 4;
imf = sprintf('/groups/jain/home/bogovicj/projects/crackSegmentation/groundTruth/closeup/img_ds%d.tif', ds)
msf = sprintf('/groups/jain/home/bogovicj/projects/crackSegmentation/groundTruth/closeup/labels_interp_smooth_ds%d.tif', ds)


i_info = imfinfo(imf);
m_info = imfinfo(msf);

N = numel(i_info);
if( N ~= numel( m_info ) )
    error('Mask and Image are different sizes');
    %    return;
end

tmp = imread(imf);
I = zeros( [ size(tmp) N ]);
M = zeros( [ size(tmp) N ]);
for i=1:N
    I(:,:,i) = imread( imf, i, 'Info', i_info );
    M(:,:,i) = imread( msf, i, 'Info', m_info );
end

clear tmp;


%% grab boundary points

% bnd = boundaryMap( (M > 0.5));
% bndi = find( bnd > 0 );

%% deal with edgels read from a file

% edgelFn = '/groups/jain/home/bogovicj/projects/crackPatching/closeup/matches_ds4/edgelMatchesNear_186.csv';
edgelFn = '/groups/jain/home/bogovicj/projects/crackPatching/closeup/matches_ds4/edgelMatches_Rad-25.0_Near_186_filt.csv'

testPt = [29.68, 108.97, 64.99];

edgelDat = csvread(edgelFn);
edgelPts = edgelDat(:, 1:3);

% be stupid for now
edgelPts = round( edgelPts );
bndi = sub2ind( size(I), edgelPts(:,1), edgelPts(:,2), edgelPts(:,3));

%% deal with all points
N = length( bndi );

% del = 5;
del = 1;

% dimrng = 1:3;
dimrng = 3;


features = cell(N,1);
hasResults = false( N, length(dimrng) );

for i = 1:N
    
    if ( mod( i, 50) == 0 )
       fprintf('working on point %d of %d\n', i, N); 
    end
    
    [x,y,z] = ind2sub( size(M), bndi(i));
    pt = [x y z];
    
    for di = 1:length(dimrng)
        d = dimrng(di);
        
        slc = getSlice( I, d, pt(d))./255;
        [feats, pts] = extractFeatures(slc, pt([1 2 3] ~= d), 'Method', 'FREAK');
        
        if( feats.NumFeatures > 0 )
            hasResults( i, di ) = true;
            features{i} = feats;
        end
        
    end
    
end

sz = size( I );
clear I M slc;

%%

fn = fullfile( SAVEPATH, 'Freak_features.mat');
save( fn );
 