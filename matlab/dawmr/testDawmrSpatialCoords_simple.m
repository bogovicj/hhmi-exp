% testDawmrSpatialCoords_simple
% 

%%

basedir = '/groups/saalfeld/home/bogovicj/projects/aniso/toys/simple';

sz = [21 21 21];
v21 = zeros( sz );
v21(11,11,11) = 1;

m21 = ones( sz, 'uint8' );

data_v21_fn = sprintf('%s/v21.h5',   basedir);
m21_fn       = sprintf('%s/m21.h5',   basedir);

if( false )
    writeH5( data_v21_fn,   v21 );
    writeH5( m21_fn,   m21, 'uint8' );
end

% centroids = eye( 125 );
% centroids = rand( 20, 125, 'single');

%% 
c_thresh_pol_kmeans = 3;
c_ave_pooling       = 0;
c_max_pooling       = 1;

patch_dim           = 5;

pooling_radius      = 0
pooling_type        = c_max_pooling;

num_clusters        = 50
num_patches_kmeans  = 1000;
num_train_instances = 0.01
num_test_instances  = Inf;

% feature_normalization = 3;
feature_normalization = 0;

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

%% v11

ds = dawmr_set( { data_v21_fn,   '', data_v21_fn   }, ...
                { m21_fn, '', m21_fn }, ...
                { m21_fn,   '', m21_fn   }, ...
                data_v21_fn );
ds.affinity_edges = [];

%%

dc = dawmr_clustering(patch_dim, num_clusters);
dp_cen  = dawmr_pooling(pooling_type, 0, [0 0 0]);
dc.add_dp(dp_cen);

%%
% dm = dawmr(ds, feature_normalization, ec_mlp_fs(mlp_init), ...
%         'v11_scriptName');

% structDm = load('/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0009_crackSegDawmr_test/exp0009_crackSegDawmr_test_results.mat', ...
%      'dm');
structDm = load('/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0008_crackSegDawmr_test/exp0008_crackSegDawmr_test_results.mat','dm'); 
dm = structDm.dm;
dm.ds = ds;

dm.sn3_mn  = {};
dm.sn3_std = {};
dm.svm_normalization = 0;

%%
% dm = dawmr_add_multilayer(dm, dc, [1]);
dm.dds.dcs.centroids = centroids;
dm.dds.dcs.num_clusters = size(centroids,1);

%%

x_test = 11;
y_test = 11;
z_test = 11;

dm.ds = ds;
dm.svm_normalization = 0;

cube_size = 1;


xs = [ x_test x_test ];
ys = [ y_test y_test ];
zs = [ z_test z_test ];

flag = 3;
data_size = dm.ds.get_data_size(flag,1);

feats = dawmr_feat_comp( dm, xs, ys, zs, data_size, flag, cube_size );

size( feats ) 
