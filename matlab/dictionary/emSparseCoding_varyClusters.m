% emSparseCoding_varyClusters
%
% dbstop if error; run_script('emSparseCoding_varyClusters','test emSparseCoding');
% dbstop if error; run_script('emSparseCoding_varyClusters','vary cluster sizes');

global SAVEPATH
global SAVEPREFIX 

%%
patch_size   = [ 9 9 9 ];
patch_size_test = [ 9 9 9 ];

num_training_patches = 10000;
num_testing_patches  = 10000;
dict_iters = 50;

% resamp_factor = [ 1 3 ];

%% data

destdir = '/groups/saalfeld/home/bogovicj/projects/aniso/2dTo3dDict/exps/recon';

training = 18;
test = 4;
ds = ds_medulla(training, test);

data_fn_train = ds.data_fn{1};
data_fn_test  = ds.data_fn{3};

%% load training patches

num_patches_tot = num_training_patches + num_testing_patches;
train_im = read_image_stack( data_fn_train );
[X, coords] = grabPatchesSimple( train_im, patch_size, num_patches_tot );
X = double(X');

X_trn = X( :, 1 : num_training_patches );
X_tst = X( :, num_training_patches+1 : end );

clear X train_im;

%%  learning params

param.iter = dict_iters;    % number of iterations

param.lambda=0.15;
param.numThreads = 1; % number of threads
param.batchsize=400;
param.verbose=false;

    
%%
num_clusters_set = [50 100 250 500 1000 1500 2000]

D_list = {};

R_list = [];
R_list_tst = [];

alpha_list = {};
alpha_list_tst = {};

i = 1;
for num_clusters = num_clusters_set
    
    param.K    = num_clusters;  % dictionary size
    
    
    %% using inria code
    % training patches should be (num_variables x num_observations)
    
    tic;
    D = mexTrainDL(X_trn,param);
    t=toc;
    fprintf('time of computation for Dictionary Learning: %f\n',t);
    
    D_list{i} = D;
    
    fprintf('Evaluating cost function for training data...\n');
    alpha_trn=mexLasso(X_trn,D,param);
    R_trn = mean( 0.5*sum(( X_trn - D * alpha_trn ).^2) ...
            + param.lambda * sum( abs(alpha_trn) ));
    
    
    R_list(i) = R_trn;
    alpha_list{i} = alpha_trn;
    
    
    fprintf('Evaluating cost function for testing data...\n');
    alpha_tst = mexLasso(X_tst,D,param);
    R_tst = mean( 0.5*sum(( X_tst - D * alpha_tst ).^2) ...
            + param.lambda * sum( abs(alpha_tst) ));
    
    R_list_tst(i) = R_tst;
    alpha_list_tst{i} = alpha_tst; 
    
    i = i + 1;
    
end

clear X_trn X_tst D alpha_trn alpha_tst;

save( sprintf( '%s/%s_results', SAVEPATH, SAVEPREFIX));
% datestr(now,30)