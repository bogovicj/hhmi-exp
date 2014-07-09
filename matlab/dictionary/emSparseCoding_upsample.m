% emSparseCoding_varyClusters
%
% dbstop if error; run_script('emSparseCoding_varyClusters','test emSparseCoding');
% dbstop if error; run_script('emSparseCoding_varyClusters','vary cluster sizes');
% dbstop if error; run_script('emSparseCoding_varyClusters','vary cluster sizes, 100 examples');
%
% ./qsubMatlabExp 'emSparseCoding_varyClusters' 'test emSparseCoding' 2
% ./qsubMatlabExp 'emSparseCoding_varyClusters' 'vary cluster sizes, 20000 examples' 4
% ./qsubMatlabExp 'emSparseCoding_varyClusters' '4k clusters, 20000 examples' 4
%
% ./qsubMatlabExp 'emSparseCoding_varyClusters' '1k clusters, 15x15x3 aniso, 10000 ex' 4

global SAVEPATH
global SAVEPREFIX 

%%
patch_size   = [ 25 25 25 ];
patch_size_test = [ 25 25 5 ];

% patch_size   = [ 15 15 3 ];
% patch_size_test = [ 15 15 3 ];

num_training_patches = 1000;
num_testing_patches  = 1000;
dict_iters = 50;

%%

dict_fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0062_emSparseCoding_varyPatchSize/exp0062_results.mat';
D_list_struct = load( dict_fn, 'D_list' );
D = D_list_struct.D_list{2};

%%  learning params

num_clusters = 50;

param.K    = num_clusters;  % dictionary size
param.iter = dict_iters;    % number of iterations

param.lambda=0.15;
param.numThreads = 4; % number of threads
param.batchsize=400;
param.verbose=false;

fprintf('building a dictionary with %d clusters\n', num_clusters);

%% data

factor = [1 1 5];
trnOnDs = 0;  % do training on downsampled volume

datdir = '/groups/saalfeld/home/bogovicj/projects/aniso/downsamp/medulla/dsdat';

training = 18;
test = 4;
ds = ds_medulla(training, test);

data_fn_train = ds.data_fn{1};
data_fn_test  = ds.data_fn{3};

mask_fn_train = ds.mask_fn{1};
mask_fn_test  = ds.mask_fn{3};

labels_fn_train = ds.labels_fn{1};
labels_fn_test  = ds.labels_fn{3};

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

% train on isotropic
if ( trnOnDs )
    ds = dawmr_set( { trnVolDsFn, '', tstVolDsFn   }, ...
                    { trnLabDsFn, '', tstLabDsFn }, ...
                    { trnMskDsFn, '', tstMskDsFn   }, ...
                    trnVolDsFn );
end

data_fn_train = ds.data_fn{1};
data_fn_test  = ds.data_fn{3};

%% load training patches

if( isempty( dict_fn ))
    
    train_im = read_image_stack( data_fn_train );
    [X_trn] = grabPatchesSimple( train_im, patch_size, num_training_patches );
    X_trn = double(X_trn');
    clear train_im;
    
    % test_im = read_image_stack( data_fn_test );
    % [X, coords] = grabPatchesSimple( test_im, patch_size_test, num_testing_patches );
    % X = double(X');
    % X_tst = X( :, num_training_patches+1 : end );
    %
    % clear X test_im
    
    % dstest = downsamplePatches( X_trn, patch_size, factor );


%% using inria code
% training patches should be (num_variables x num_observations)


    tic;
    D = mexTrainDL(X_trn,param);
    t=toc;
    fprintf('time of computation for Dictionary Learning: %f\n',t);
end

% % SKIP EVALUATION
% fprintf('Evaluating cost function for training data...\n');
% alpha_trn=mexLasso(X_trn,D,param);
% R_trn = mean( 0.5*sum(( X_trn - D * alpha_trn ).^2) ...
%     + param.lambda * sum( abs(alpha_trn) ));

% fprintf('Evaluating cost function for testing data...\n');
% alpha_tst = mexLasso(X_tst,D,param);
% R_tst = mean( 0.5*sum(( X_tst - D * alpha_tst ).^2) ...
%     + param.lambda * sum( abs(alpha_tst) ));


% i = i + 1;

%% downsample dict

D_ds = downsamplePatches( D, patch_size, factor );

%% test reconstruction on test data

test_im = read_image_stack( data_fn_test );
[X_tst] = grabPatchesSimple( test_im, patch_size, num_testing_patches );
X_tst = double(X_tst');
clear test_im;

%%
fprintf('Downsampling test patches...\n');
X_tst_ds = downsamplePatches( X_tst, patch_size, factor );

fprintf('Computing new representation for aniso data...\n');
alpha_tst = mexLasso( X_tst_ds, D_ds, param);

X_tst_ds_us = D * alpha_tst;

%%

i = 1;
orig_patch = reshape( X_tst(:,1), patch_size );
re_patch = reshape( X_tst_ds_us(:,1), patch_size );

BrowseComponents('iii', orig_patch, re_patch, orig_patch - re_patch);


%% save
save( sprintf( '%s/%s_results', SAVEPATH, SAVEPREFIX));
% datestr(now,30)