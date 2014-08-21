% genFromTif

%% set up
srcdir = '/nobackup/saalfeld/john/forMariaLuisa';
dstdir = '/groups/saalfeld/home/bogovicj/projects/davi/mariaAnnot';

%% params

do_full = 1;
id = '1';
prefix = '140718_elastic_registered_AL_Maria-Luisa_ground_truth_v04';

confThresholds = [0.55 0.75 1.0 ]; % lo mid hi
confLabels      = {'lo','med','hi'};
nConfs = length( confLabels );

norm_by_slice = 1;
pad_amt       = 4;

%%
lab_list = { 'neuron', ...
             'mito',   ...
             'glia' };
         
%% preliminaries

nLabels = length( lab_list );

%% image to5

imrgb = readMultiTiff( fullfile( srcdir, id, 'img_rgb.tif'));
sz = size(imrgb);
sz = sz(1:3)
% mx = max(imrgb(:));

% generate mask
msk = ~((imrgb(:,:,:,1)==0) & (imrgb(:,:,:,2)==255) & (imrgb(:,:,:,3)==0));

% normalize
if( norm_by_slice )
    im = imrgb(:,:,:,1);
    for z = 1:size(imrgb,3)
        slc = im(:,:,z);
        slcm = msk(:,:,z);
        
        mn = mean( slc(slcm) );
        sd = std( slc(slcm) );
        
        im(:,:,z) = (slc - mn) ./ sd;
    end
else
    im = imrgb(:,:,:,1);
    
    mn = mean( im(msk) );
    sd = std ( im(msk) );
    im = ( im - mn )./ sd;
end

clear imrgb slc slcm;

%% label to h5

lbvol = zeros( [ sz nLabels ] );

for i = 1:nLabels
    for c = 1:nConfs
        thisim =  readMultiTiff(  fullfile( srcdir, id, [lab_list{i},'_', confLabels{c}, '.tif']) );
        unique( thisim(:))
        lbvol(:,:,:,i) = max( lbvol(:,:,:,i), confThresholds(c) .* thisim);
    end
end
clear thisim;   

%% update mask

msk = repmat( msk, [1 1 1 nLabels] );

%% pad

if( ~isempty( pad_amt))
    im    = pad_volume( im, 4, 0 );
    lbvol = pad_volume( lbvol, 4, 0 );
    msk   = pad_volume( msk, 4, 0 );
end

%% write everything

system(sprintf('mkdir -p %s',fullfile(dstdir,id)))

writeH5( fullfile(dstdir,id, sprintf('%s_im.h5', prefix)),  im );
writeH5( fullfile(dstdir,id, sprintf('%s_lb.h5', prefix)),  lbvol );
writeH5( fullfile(dstdir,id, sprintf('%s_msk.h5', prefix)), uint8(msk), 'uint8');

