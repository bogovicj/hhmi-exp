%% dictionaryLearningForAnsio
% run_script('dictionaryLearningForAnsio', 'patchR 7, clusters 50');
% run_script('dictionaryLearningForAnsio', 'iso - patchR 5, clusters 1000, scales [1]');
% run_script('dictionaryLearningForAnsio', 'aniso z2 - patchR 5 5 3, clusters 1000, scales [1]');

global SAVEPATH
global SAVEPREFIX

script_name = mfilename;

%%
factor = [1 1 2];
trnOnDs = 1;  % do training on downsampled volume

datdir = '/groups/saalfeld/home/bogovicj/projects/aniso/downsamp/medulla/dsdat';

trnVolumeFn = ['/groups/jain/home/jainv/datasets/medulla_oct12/' ...
                'validation/im_normalized_0mean.h5'];

tstVolumeFn = ['/groups/jain/home/jainv/datasets/medulla_oct12/' ...
                'testing/im_normalized_0mean.h5'];
    
%
training = 18;
base_dir_ira = sprintf('/groups/jain/home/huangg/research/exp/medulla_ira_training%d/', training);
labels_fn_train = [base_dir_ira 'labels.h5'];
mask_fn_train   = [base_dir_ira 'mask.h5'];

test = 3;
base_dir_gbh = sprintf(['/groups/jain/home/huangg/research/exp/medulla_ira_validation%d/'], test);
labels_fn_test = [base_dir_gbh 'labels.h5'];
mask_fn_test   = [base_dir_gbh 'mask.h5'];

[~,trnVolName] = fileparts( trnVolumeFn );
[~,tstVolName] = fileparts( tstVolumeFn );
trnVolDsFn = fullfile(datdir, sprintf('trn_%s_ds%d-%d-%d.h5',trnVolName,factor));
tstVolDsFn = fullfile(datdir, sprintf('tst_%s_ds%d-%d-%d.h5',tstVolName,factor));

[~,trnMskName] = fileparts( mask_fn_train );
[~,tstMskName] = fileparts( mask_fn_test );
trnMskDsFn = fullfile(datdir, sprintf('trn_%s_ds%d-%d-%d.h5',trnMskName,factor));
tstMskDsFn = fullfile(datdir, sprintf('tst_%s_ds%d-%d-%d.h5',tstMskName,factor));

[~,trnLabName] = fileparts( labels_fn_train );
[~,tstLabName] = fileparts( labels_fn_test );
trnLabDsFn = fullfile(datdir, sprintf('trn_%s_ds%d-%d-%d.h5',trnLabName,factor));
tstLabDsFn = fullfile(datdir, sprintf('tst_%s_ds%d-%d-%d.h5',tstLabName,factor));

if( ~exist(trnVolDsFn,'file') )
    fprintf('writing test volume %d-%d-%d\n',factor);
    success = downsampWriteH5( trnVolumeFn, trnVolDsFn, factor );
end
if( ~exist(trnMskDsFn,'file') )
    fprintf('writing training mask %d-%d-%d\n',factor);
    success = downsampWriteH5( mask_fn_train, trnMskDsFn, factor, 'or' );
end
if( ~exist(trnLabDsFn,'file') )
    fprintf('writing training labels %d-%d-%d\n',factor);
    success = downsampWriteH5( labels_fn_train, trnLabDsFn, factor, 'avg' );
end

if( ~exist(tstVolDsFn,'file') )
    fprintf('writing test volume %d-%d-%d\n',factor);
    success = downsampWriteH5( tstVolumeFn, tstVolDsFn, factor );
end
if( ~exist(tstMskDsFn,'file') )
    fprintf('writing test mask %d-%d-%d\n',factor);
    success = downsampWriteH5( mask_fn_test, tstMskDsFn, factor, 'or' );
end
if( ~exist(tstVolDsFn,'file') )
    fprintf('writing test labels %d-%d-%d\n',factor);
    success = downsampWriteH5( labels_fn_test, tstLabName, factor, 'avg' );
end


% train on isotropic
if ( trnOnDs )
    ds = dawmr_set( { trnVolDsFn,   '', ''   }, ...
                    { trnLabDsFn, '', '' }, ...
                    { trnMskDsFn,   '', ''   }, ...
                    trnVolDsFn );
else
    ds = dawmr_set( { trnVolumeFn,   '', ''   }, ...
                    { labels_fn_train, '', '' }, ...
                    { mask_fn_train,   '', ''   }, ...
                    trnVolumeFn );
end
      
ds.affinity_edges = [];

%% specify unsupervised architecture
c_thresh_pol_kmeans = 3;
c_ave_pooling       = 0;
c_max_pooling       = 1;
patch_dim           = [5 5 3]

pooling_radius      = 2;
pooling_type        = c_max_pooling;

num_clusters        = 1000
num_patches_kmeans  = 10000;
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
dm = dawmr_add_multilayer( dm, dc, [1]);
% dm = dawmr_add_multilayer( dm, dc, [1 2]);
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
