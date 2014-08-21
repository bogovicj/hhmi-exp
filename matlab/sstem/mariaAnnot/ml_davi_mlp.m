% quick start sample script
% maria luisa's annotations on davis data
%
% dbstop if error; run_script('ml_davi_mlp','dictionary only, 1k clusters');


global SAVEPATH
global SAVEPREFIX

%% exp params

saved_dict_fn = '';

do_classifier   = 1;

do_foveated = 1;
do_cn       = 0;

num_workers = 50;

%% define the data set

training = 1;
test = [];
ds = ds_medulla(training, test);

%% specify unsupervised architecture

alphas    = 0;

c_thresh_pol_kmeans = 3;
c_ave_pooling       = 0;
c_max_pooling       = 1;
patch_dim           = [7 7 3];

encoding_type      = c_thresh_pol_kmeans;
pooling_radius     = [2 2 1];
pooling_type       = c_max_pooling;

num_clusters        = 1000;
num_patches_kmeans  = 10000;
num_train_instances = Inf;
num_test_instances  = Inf;

feature_normalization = 3;
% feature_normalization = 0;

if (do_foveated)
    dc = dawmr_clustering(patch_dim, num_clusters, do_cn, ...
        [], [],[],[], [],[], 0);
    dc = dc_foveated(dc, pooling_radius, pooling_type, encoding_type );
else
    dc = dawmr_clustering(patch_dim, num_clusters);
    dp_cen  = dawmr_pooling(pooling_type, 0, [0 0 0]);
    dc.add_dp(dp_cen);
end

downsampling_type = 1;


%% specify mlp parameters

mlp_init = mlp( 100 );
mlp_init.num_updates_default = 5e5;
mlp_init.minibatch_size_default = 40;
mlp_init.use_gpu = 1;
mlp_init.pos_ratio = 0.5;
mlp_init.margin = 0.4;
mlp_init.eta_w_start = [0.02 0.02];
mlp_init.eta_b_start = [0.02 0.02];
mlp_init.eta_b_start = [0.02 0.02];
mlp_init.loss_func = 1;

% dm.end_classifier = ec_mlp_fs(mlp_init);

%%

if( isempty(saved_dict_fn) )
    
    dm = dawmr(ds, 3, ec_mlp_fs(mlp_init), script_name);
    dm = dawmr_add_multilayer(dm, dc, [1]);
    dm.learn_features(num_patches_kmeans, 2, [],[],[], 50);
    
    if(~isempty(SAVEPATH))
        save(sprintf('%s/%s_%s_learnedFeatures.mat', ...
            SAVEPATH, SAVEPREFIX, script_name), ...
            'dm','-v7.3');
    end
end


if( do_classifier )
    
    ts = tic;
    
    % train 
    [accs_train, labels_gt_train, ~, labels_pd_train] = ...
        dm.classifier(1, num_train_instances, ...
            [],[],[],[], num_workers);
      
    
    if( ~isempty(test))
        % test 
        
        [accs, labels_gt, ~, labels_pd, aucs] = ...
            dm.classifier(3, num_test_instances, ...
            [],[],[],[], num_workers);
    end
    
    te = toc(ts);
    
    fprintf('first model\n');
    dawmr_set_print_stats(accs_train, accs, aucs, te, [], ...
        labels_gt_train, labels_pd_train, ...
        labels_gt, labels_pd, ...
        SAVEPREFIX, script_name, dm, num_clusters);
    
    % results save
    if(~isempty(SAVEPATH))
        save(sprintf('%s/%s_%s_finalResults.mat', ...
            SAVEPATH, SAVEPREFIX, script_name), ...
            '-v7.3');
    end
    
end
