%% dictionaryLearningForAnsio
% run_script('dictionaryLearningForAnsio', 'patchR 7, clusters 50');
% run_script('dictionaryLearningForAnsio', 'iso - patchR 5, clusters 1000, scales [1]');
% run_script('dictionaryLearningForAnsio', 'aniso z2 - patchR 5 5 3, clusters 1000, scales [1]');
%
% run_script('dictionaryLearningForAnsio', 'MEDULLA SAMPLE, iso - patch 5x5x5, clusters 1000, scales [1]');

global SAVEPATH
global SAVEPREFIX

script_name = mfilename;

%%
factor = [1 1 2];
trnOnDs = 1;  % do training on downsampled volume

% datdir = '/groups/saalfeld/home/bogovicj/projects/aniso/downsamp/sample_medulla/dsdat';
% base_dir = '/groups/saalfeld/home/bogovicj/dev/dawmr/dawmr_lib_public/projects/sample_medulla';
% 
% data_fn_train   = sprintf('%s/medulla_sub1_data.h5',   base_dir);
% mask_fn_train   = sprintf('%s/medulla_sub1_mask_1.h5',   base_dir);
% labels_fn_train = sprintf('%s/medulla_sub1_labels_1.h5', base_dir);
% 
% data_fn_test    = sprintf('%s/medulla_sub2_data.h5',   base_dir);
% mask_fn_test    = sprintf('%s/medulla_sub2_mask_1.h5',   base_dir);
% labels_fn_test  = sprintf('%s/medulla_sub2_labels_1.h5', base_dir);
   
training = 18;
test = 4;
ds = ds_medulla(training, test);

data_fn_train = ds.data_fn{1};
data_fn_test  = ds.data_fn{3};

mask_fn_train = ds.mask_fn{1};
mask_fn_test  = ds.mask_fn{3};

labels_fn_train = ds.labels_fn{1};
labels_fn_test  = ds.labels_fn{3};

%%

[~,trnVolName] = fileparts( data_fn_train );
[~,tstVolName] = fileparts( data_fn_test );
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
    fprintf('writing training volume %d-%d-%d\n',factor);
    success = downsampWriteH5( data_fn_train, trnVolDsFn, factor );
end
if( ~exist(trnMskDsFn,'file') )
    fprintf('writing training mask %d-%d-%d\n',factor);
    success = downsampWriteH5( mask_fn_train, trnMskDsFn, factor );
end
if( ~exist(trnLabDsFn,'file') )
    fprintf('writing training labels %d-%d-%d\n',factor);
    success = downsampWriteH5( labels_fn_train, trnLabDsFn, factor );
end

if( ~exist(tstVolDsFn,'file') )
    fprintf('writing test volume %d-%d-%d\n',factor);
    success = downsampWriteH5( data_fn_test, tstVolDsFn, factor );
end
if( ~exist(tstMskDsFn,'file') )
    fprintf('writing test mask %d-%d-%d\n',factor);
    success = downsampWriteH5( mask_fn_test, tstMskDsFn, factor );
end
if( ~exist(tstLabDsFn,'file') )
    fprintf('writing test labels %d-%d-%d\n',factor);
    success = downsampWriteH5( labels_fn_test, tstLabDsFn, factor );
end


% train on isotropic
if ( trnOnDs )
    ds = dawmr_set( { trnVolDsFn,   '', tstVolDsFn   }, ...
                    { trnLabDsFn, '', tstLabDsFn }, ...
                    { trnMskDsFn,   '', tstMskDsFn   }, ...
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
patch_dim           = [5 5 5]

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
