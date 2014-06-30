%% upsampPatchDawmrCentroids

% load '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0015_dictionaryLearningForAnsio/exp0015_dictionaryLearningForAnsio_results.mat';
load /groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0019_dictionaryLearningForAnsio/exp0019_dictionaryLearningForAnsio_results.mat

%%
cluster_cents = dm.dds.dcs(1).centroids;

rng('default');
rng(1);

N = 20;
sample_pts = randi( 25, [ N 3 ] );

factor = [1 1 5];


for i = 1 : size(cluster_cents,1)
    cents_rs = reshape( cluster_cents(i,:), dc.patch_dim );
%     cents_rs_sm = resampleRatioCenter( cents_rs, r );
    cents_rs_sm_il = downSampleGaussianImglib( toImglib(cents_rs), factor );
    cents_rs_sm = net.imglib2.util.ImgOps.toFloatArray3d( cents_rs_sm_il ); 
    clear cents_rs_sm_il
    
    if ( i == 1 )
        patch_size_new = size( cents_rs_sm );
        cents_sm = zeros( size( cluster_cents,1), numel(cents_rs_sm));
    end
    cents_sm( i, : ) = cents_rs_sm(:)';
end

