% crackSegDawmr_test
% 
% run_script('crackSegDawmr_test', 'qs medulla test');
% run_script('crackSegDawmr_test', 'ds medulla test, feat norm 3');
%
% run_script('crackSegDawmr_test', 'ds_medulla: downsample 1x1x2 ds dist, ds dat');
% run_script('crackSegDawmr_test', 'ds_medulla: downsample 1x1x2 patch_5x5x3 ds dist, ds dat');
% run_script('crackSegDawmr_test', 'ds_medulla: patch_5x5x5 learn dictionary only');
% run_script('crackSegDawmr_test', 'ds_medulla: patch_5x5x5, iso dictionary, iso data');
%
% run_script('crackSegDawmr_test', 'ds_medulla: patch_5x5x5, iso dictionary downsampled, aniso data');
% run_script('crackSegDawmr_test', 'debug dictionary downsampling');
%
% run_script('crackSegDawmr_test', 'a test of crack segmentation 2');
% run_script('crackSegDawmr_test', 'crack segmentation legit 200k iters, 200 hu');
% run_script('crackSegDawmr_test', 'crack segmentation legit 20k iters, 50 clusters, 50 hu, no norm');
% run_script('crackSegDawmr_test', 'test of aniso data');

global SAVEPATH
global SAVEPREFIX

script_name = mfilename;

%% distributed computing parameters
num_workers = 50;

trnOnDs = 1;
factor = [ 1 1 2 ];

dict_exp = 'exp0035_crackSegDawmr_test'; % experiment containing the dictionary to use

dm_dict = [];
if( ~isempty(dict_exp) )
    dictFn = sprintf( '%s/../%s/%s_results.mat', SAVEPATH, dict_exp, dict_exp );
    tmp = load( dictFn, 'dm');
    dm_dict = tmp.dm;
    clear tmp;
end

%% define the data set
% % base_dir = '/groups/saalfeld/home/bogovicj/projects/crackSegmentation/groundTruth';
% base_dir = sprintf('%s/projects/sample_medulla', DAWMRLIBPATH);
% % 
% % data_fn_train   = sprintf('%s/closeup/closeup_img_ds4.h5',   base_dir);
% % mask_fn_train   = sprintf('%s/closeup/closeup_msk_ds4.h5',   base_dir);
% % labels_fn_train = sprintf('%s/closeup/closeup_lbls_ds4.h5', base_dir);
% % 
% data_fn_train   = sprintf('%s/medulla_sub1_data.h5',   base_dir);
% mask_fn_train   = sprintf('%s/medulla_sub1_mask_1.h5',   base_dir);
% labels_fn_train = sprintf('%s/medulla_sub1_labels_1.h5', base_dir);
% % 
% % data_fn_test    = sprintf('%s/closeup2/closeup2_img_ds4.h5',   base_dir);
% % mask_fn_test    = sprintf('%s/closeup2/closeup2_msk_ds4.h5',   base_dir);
% % labels_fn_test  = sprintf('%s/closeup2/closeup2_emptyLabels_ds4.h5', base_dir);
% % 
% data_fn_test    = sprintf('%s/medulla_sub2_data.h5',   base_dir);
% mask_fn_test    = sprintf('%s/medulla_sub2_mask_1.h5',   base_dir);
% labels_fn_test  = sprintf('%s/medulla_sub2_labels_1.h5', base_dir);
% 
% ds = dawmr_set( { data_fn_train,   '', data_fn_test   }, ...
%                 { labels_fn_train, '', labels_fn_test }, ...
%                 { mask_fn_train,   '', mask_fn_test   }, ...
%                 data_fn_train );
% ds.affinity_edges = [];

%%
training = 18;
test = 4;
ds_tmp = ds_medulla(training, test) 
% ds.affinity_edges = [];

data_fn_train = ds_tmp.data_fn{1};
data_fn_test  = ds_tmp.data_fn{3};

mask_fn_train = ds_tmp.mask_fn{1};
mask_fn_test  = ds_tmp.mask_fn{3};

labels_fn_train = ds_tmp.labels_fn{1};
labels_fn_test  = ds_tmp.labels_fn{3};

%%
if ( trnOnDs )
    
    datdir = '/groups/saalfeld/home/bogovicj/projects/aniso/downsamp/medulla/dsdat';
    
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
    
    ds = dawmr_set( { trnVolDsFn,   '', tstVolDsFn   }, ...
        { trnLabDsFn, '', tstLabDsFn }, ...
        { trnMskDsFn,   '', tstMskDsFn   }, ...
        trnVolDsFn );
    
    data_fn_train = ds.data_fn{1};
    data_fn_test  = ds.data_fn{3};
    
    mask_fn_train = ds.mask_fn{1};
    mask_fn_test  = ds.mask_fn{3};
    
    labels_fn_train = ds.labels_fn{1};
    labels_fn_test  = ds.labels_fn{3};
    
else
    ds = ds_tmp;
end

ds

%% specify unsupervised architecture
c_thresh_pol_kmeans = 3;
c_ave_pooling       = 0;
c_max_pooling       = 1;

patch_dim           = [5 5 5]

pooling_radius      = 2
pooling_type        = c_max_pooling;

num_clusters        = 50
num_patches_kmeans  = 1000;
num_train_instances = 0.01
num_test_instances  = Inf;

feature_normalization = 3;
% feature_normalization = 0;

dc = dawmr_clustering(patch_dim, num_clusters);
% dc = dc_foveated(dc, pooling_radius, pooling_type, ...
%                  c_thresh_pol_kmeans);

dp_cen  = dawmr_pooling(pooling_type, 0, [0 0 0]);
dc.add_dp(dp_cen);

%% specify mlp parameters
mlp_init = mlp(100);
% mlp_init.num_updates_default = 2e5;
mlp_init.num_updates_default = 2e4;
mlp_init.minibatch_size_default = 40;
mlp_init.use_gpu = 1;
mlp_init.pos_ratio = 0.5;
mlp_init.margin = 0.4;
mlp_init.eta_w_start = [0.02 0.02];
mlp_init.eta_b_start = [0.02 0.02];
mlp_init.loss_func = 1;

%% segmentation parameters
thresh           = [0.5 0.7 0.8 0.9 0.95 0.99];
watershed_thresh = [0];


%% set up model, do learning
ts = tic;
if( isempty(dm_dict))
    dm = dawmr(ds, feature_normalization, ec_mlp_fs(mlp_init), script_name);
    dm = dawmr_add_multilayer(dm, dc, [1]);
    
    fprintf('learning features\n');
    dm.learn_features(num_patches_kmeans, 2, [],[],[],50);
    fprintf('done\n');

else
    fprintf('using pre-computed features\n');
    dm = resampleDawmrClusters( dm_dict, factor );
end

[accs_train, labels_gt_train, ~, labels_pd_train] = ...
    dm.classifier(1, num_train_instances, ...
                  [],[],[],[], num_workers);
[accs, labels_gt, ~, labels_pd, aucs] = ...
    dm.classifier(3, num_test_instances, ...
                  [],[],[],[], num_workers);
te = toc(ts);

fprintf('first model\n');
dawmr_set_print_stats(accs_train, accs, aucs, te, [], ...
                      labels_gt_train, labels_pd_train, ...
                      labels_gt, labels_pd, ...
                      SAVEPREFIX, script_name, dm, num_clusters);

% % temporary save
if(~isempty(SAVEPATH))
  save(sprintf('%s/%s_%s_results.mat', ...
               SAVEPATH, SAVEPREFIX, script_name), ...
       '-v7.3');
end

% do full stack inference
% dm.infer([],[],...
%          sprintf('%s/%s_%s_sub2_infer.h5', SAVEPATH, SAVEPREFIX, script_name), ...
%          [],[],[],[20 20 20],[], 3);
% dm.infer([],[],...
%          sprintf('%s/%s_%s_sub1_infer.h5', SAVEPATH, SAVEPREFIX, script_name), ...
%          [],[],[],[20 20 20],[], 1);
