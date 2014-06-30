%% dictionaryLearningForAnsio
% run_script('dictionaryLearningForAnsio', 'patchR 7, clusters 50');
% run_script('dictionaryLearningForAnsio', 'patchR 5, clusters 100, ms [1 2]');

global SAVEPATH
global SAVEPREFIX

script_name = mfilename;

%%
factor = [1 1 5];

destdir = '/groups/saalfeld/home/bogovicj/projects/aniso/downsamp/medulla';

trnVolumeFn = ['/groups/jain/home/jainv/datasets/medulla_oct12/' ...
                'validation/im_normalized_0mean.h5'];

tstVolumeFn = ['/groups/jain/home/jainv/datasets/medulla_oct12/' ...
                'testing/im_normalized_0mean.h5'];
            
[~,trnVolName] = fileparts( trnVolumeFn );
[~,tstVolName] = fileparts( tstVolumeFn );

trnVolDsFn = fullfile(destdir, sprintf('trn_%s_ds%d-%d-%d.h5',trnVolName,factor));
tstVolDsFn = fullfile(destdir, sprintf('tst_%s_ds%d-%d-%d.h5',tstVolName,factor));

%
training = 18;
base_dir_ira = sprintf('/groups/jain/home/huangg/research/exp/medulla_ira_training%d/', training);
labels_fn_train = [base_dir_ira 'labels.h5'];
mask_fn_train   = [base_dir_ira 'mask.h5'];

test = 3;
base_dir_gbh = sprintf(['/groups/jain/home/huangg/research/exp/medulla_ira_validation%d/'], test);
labels_fn_test = [base_dir_gbh 'labels.h5'];
mask_fn_test   = [base_dir_gbh 'mask.h5'];

ds = dawmr_set( { trnVolumeFn,   '', tstVolumeFn   }, ...
                { labels_fn_train, '', labels_fn_test }, ...
                { mask_fn_train,   '', mask_fn_test   }, ...
                trnVolumeFn );
ds.affinity_edges = [];

%% downsample z

% trnVol = read_image_stack( trnVolumeFn );
% trnvolil = toImglib( trnVol ); 
% clear trnVol
% trnVolilds = downSampleGaussianImglib( trnvolil, factor ); 
% clear trnvolil
% trnVolDs = net.imglib2.util.ImgOps.toFloatArray3d( trnVolilds ); 
% clear trnVolilds
% writeH5( trnVolDsFn, trnVolDs); 
% clear trnVolDs
% 
% tstVol = read_image_stack( tstVolumeFn );
% tstvolil = toImglib( tstVol );
% clear tstVol 
% tstVolilds = downSampleGaussianImglib( tstvolil, factor );
% clear tstvolil
% tstVolDs = net.imglib2.util.ImgOps.toFloatArray3d( tstVolilds );
% clear tstVolilds
% writeH5( tstVolDsFn, tstVolDs ); 
% clear tstVolDs
% 
% % BrowseComponents('ii', trnVol, tstVol);


%% specify unsupervised architecture
c_thresh_pol_kmeans = 3;
c_ave_pooling       = 0;
c_max_pooling       = 1;
patch_dim           = 5

pooling_radius      = 2;
pooling_type        = c_max_pooling;

num_clusters        = 100
num_patches_kmeans  = 1000;
num_train_instances = Inf;
num_test_instances  = Inf;

% feature_normalization = 3;
feature_normalization = 0;

dc = dawmr_clustering(patch_dim, num_clusters);

dp_cen  = dawmr_pooling(pooling_type, 0, [0 0 0]);
dc.add_dp(dp_cen);

downsampling_type = 1;


mlp_init = mlp(100); % doesn't matter


%% set up model, do learning
ts = tic;
dm = dawmr(ds, feature_normalization, ec_mlp_fs(mlp_init), script_name);
dm = dawmr_add_multilayer( dm, dc, [1 2]);
% dm = dawmr_add_multilayer( dm, dc, scales );

fprintf('learning features\n');
cube_size = 200;
dm.learn_features( num_patches_kmeans, 2, [], [], [], cube_size );
fprintf('done\n');

% % save
if(~isempty(SAVEPATH))
  save(sprintf('%s/%s_%s_results.mat', ...
               SAVEPATH, SAVEPREFIX, script_name), ...
       '-v7.3');
end
