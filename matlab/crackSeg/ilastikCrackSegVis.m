%% ilastikCrackSegVis
% Visualize crack segmentation performed by ilastik random forest classifier
% Assumes that 

%% parameters
basedir = '/groups/saalfeld/home/bogovicj/projects/crackSegmentation';
traintag = 'closeup'
ilastikProj = 'closeup_ds4_2';
xvaltags = {'closeup2', 'closeup3'};

ds = 4; % the downsample factor

% gray colormap
gmap = repmat( linspace(0,1,128)', 1, 3 );

% prob colormap
cmap = diverging_map(linspace(0,1,128),[0 0 1],[1 0 0]);
 
%% set up

% the ilastik project directory
ilastikDir = fullfile( basedir, 'ilastik', ilastikProj);

% the image file name
imfn = fullfile( basedir, 'groundTruth', ...
        traintag, sprintf('%s_img_ds%d.tif', traintag, ds));
segfn = fullfile( ilastikDir, sprintf('%s_img_ds%d_export.h5',traintag, ds));

% read the segmentation
seg = h5read(segfn, '/exported_data');

% read the image
im = readMultiTiff( imfn );

%%

% so ugly!
maskSeg = squeeze(seg(2,:,:,:));
maskSeg = flipdim(permute(flipdim(permute(maskSeg, [3 2 1]), 1), [2 1 3]),2);

% BrowseComponents('ii', im, maskSeg);

%% 

slc = 72;

imSlc = im(:,:,slc)./255;
mkSlc = maskSeg(:,:,slc);

immk = overlayImagesGen( imSlc, 127.*mkSlc, [0 1.1], gmap, cmap, 0.33 );
% immk = overlayImages( imSlc, mkSlc, gmap, cmap, [0 1], 0.2, 0);
figure('color','w');
image( immk );
axis equal; axis off;

%%

[~,segname]=fileparts(segfn);
fnout = fullfile( ilastikDir, 'vis', sprintf('%s-%s_slc%d.png',ilastikProj,segname,slc));
export_fig( fnout, '-m2' );
