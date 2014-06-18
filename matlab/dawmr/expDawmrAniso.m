%% expDawmrAniso
% experiment with dawmr clustering / segmentation using aniso data

%% open an older model

savedDatFn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0008_crackSegDawmr_test/exp0008_crackSegDawmr_test_results.mat';

load( savedDatFn );

%% generate a downsampled image
basedir = '/groups/saalfeld/home/bogovicj/projects/aniso/downsamp';
[~,data_fn_train_ds] = fileparts(strrep(dm.ds.data_fn{1},'.h5','_lr2.h5'))
[~,mask_fn_train_ds] = fileparts(strrep(dm.ds.mask_fn{1},'.h5','_lr2.h5'))
[~,lbl_fn_train_ds ] = fileparts(strrep(dm.ds.labels_fn{1},'.h5','_lr2.h5'))

data_fn_train_ds = [fullfile(basedir, data_fn_train_ds) '.h5']
mask_fn_train_ds = [fullfile(basedir, mask_fn_train_ds) '.h5']
lbl_fn_train_ds  = [fullfile(basedir, lbl_fn_train_ds) '.h5']

%%

if( ~ exist( data_fn_train_ds, 'file' ) )
    disp('writing data');
    im = read_image_stack( dm.ds.data_fn{1} );
    msk = read_image_stack( dm.ds.mask_fn{1} );
    lbl = read_image_stack( dm.ds.labels_fn{1} );
    
    im_ds  = reduceEffectiveResolution( im, 2 );
    msk_ds = reduceEffectiveResolution( msk, 2 );
    lbl_ds = reduceEffectiveResolution( lbl, 2 );
    
    writeH5( data_fn_train_ds, im_ds, 'single', [10 10 10]);
    writeH5( mask_fn_train_ds, msk_ds, 'uint8', [10 10 10]);
    writeH5( lbl_fn_train_ds,  lbl_ds, 'uint8', [10 10 10]);
end

%% learn features

ds = dawmr_set( { data_fn_train,   '', data_fn_train_ds  }, ...
                { labels_fn_train, '', lbl_fn_train_ds }, ...
                { mask_fn_train,   '', mask_fn_train_ds }, ...
                data_fn_train );
ds.affinity_edges = [];

dm.ds = ds;
this = dm;

%%

fnout = sprintf('%s/infer_%s.h5', basedir, date );
dm.infer([],[],...
    fnout, ...
    [],[],[],[20 20 20],[], 3);
     

%% open the true and new results

savedRes = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0008_crackSegDawmr_test/exp0008_crackSegDawmr_test_sub1_infer.h5';
seg = read_image_stack( savedRes );
segnew = read_image_stack( fnout );

%%

BrowseComponents('iii', seg, segnew, (seg - segnew))

