%% bpFromSavedDict
% run_script('bpFromSavedDict', 'test');

global SAVEPATH
global SAVEPREFIX

script_name = mfilename;
%% distributed computing parameters
num_workers = 75;

%%

% SAVEPATH = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0025_dictionaryLearningForAnsio';

dictExp = '0026';
dictFn = sprintf( '%s/../exp%s_dictionaryLearningForAnsio/exp%s_dictionaryLearningForAnsio_results.mat', SAVEPATH, dictExp, dictExp );
exist( dictFn, 'file')

dictRes = load( dictFn );
dm = dictRes.dm;
clear dictRes

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

test = 4;
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


% train on isotropic

ds = dawmr_set( { trnVolDsFn, '', tstVolDsFn   }, ...
                { trnLabDsFn, '', tstLabDsFn }, ...
                { trnMskDsFn, '', tstMskDsFn   }, ...
                trnVolDsFn );

ds.affinity_edges = [];

%% specify supervised architecture

% feature_normalization = 3;
feature_normalization = 0;

downsampling_type = 1;

num_train_instances = Inf;
num_test_instances  = Inf;

mlp_init = mlp(100);
mlp_init.num_updates_default = 5e5;
mlp_init.minibatch_size_default = 40;
mlp_init.use_gpu = 1;
mlp_init.pos_ratio = 0.5;
mlp_init.margin = 0.4;
mlp_init.eta_w_start = [0.02 0.02];
mlp_init.eta_b_start = [0.02 0.02];
mlp_init.loss_func = 1;


%% set up model 

dm.end_classifier = ec_mlp_fs(mlp_init);

% we have the features learned here already
dm.sn3_mn  = {};
dm.sn3_std = {};
dm.svm_normalization = 0;

dm.ds = ds;

%%

ts = tic;

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

% % save the mlp
if(~isempty(SAVEPATH))
  save(sprintf('%s/%s_%s_results.mat', ...
               SAVEPATH, SAVEPREFIX, script_name), ...
       '-v7.3');
end

% dont bother with the full stack inference for now ...
