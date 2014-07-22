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

%%


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

%% label to h5

lbvol = zeros( [ sz nLabels ] );

for i = 1:nLabels
    fn =  fullfile( srcdir, id, [lab_list{i},'_conf.tif']);
    if( ~exist( fn, 'file'))
        for c = 1:nConfs
            thisim =  readMultiTiff(  fullfile( srcdir, id, [lab_list{i},'_', confLabels{c}, '.tif']) );
            unique( thisim(:))
            lbvol(:,:,:,i) = lbvol(:,:,:,i) + confThresholds(c) .* thisim;
        end
    end
end
clear thisim;   

%% generate mask

mkvol = ones( size(lbvol), 'uint8' );
mkvol(:,:,11,:) = 0; % slice 11 missing

%% write everything

system(sprintf('mkdir -p %s',fullfile(dstdir,id)))

writeH5( fullfile(dstdir,id, sprintf('%s_im.h5', prefix)),  imvol);
writeH5( fullfile(dstdir,id, sprintf('%s_lb.h5', prefix)),  lbvol, 'uint8');
writeH5( fullfile(dstdir,id, sprintf('%s_msk.h5', prefix)), mkvol, 'uint8');

