%% expDawmrAniso
% experiment with dawmr clustering / segmentation using aniso data

%% open an older model

savedDatFn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0008_crackSegDawmr_test/exp0008_crackSegDawmr_test_results.mat';
load( savedDatFn );

% for dsFactor = [ 3 5 7 9]
%     
%     ds_suffix =  ['_lr', num2str(dsFactor) ,'.h5']
%     
%     
%     % generate a downsampled image
%     basedir = '/groups/saalfeld/home/bogovicj/projects/aniso/downsamp';
%     [~,data_fn_train_ds] = fileparts(strrep(dm.ds.data_fn{1},  '.h5', ds_suffix ));
%     [~,mask_fn_train_ds] = fileparts(strrep(dm.ds.mask_fn{1},  '.h5', ds_suffix ));
%     [~,lbl_fn_train_ds ] = fileparts(strrep(dm.ds.labels_fn{1},'.h5', ds_suffix ));
%     
%     data_fn_train_ds = [fullfile(basedir, data_fn_train_ds) '.h5']
%     mask_fn_train_ds = [fullfile(basedir, mask_fn_train_ds) '.h5']
%     lbl_fn_train_ds  = [fullfile(basedir, lbl_fn_train_ds) '.h5']
%     
%     ds = dawmr_set( { data_fn_train,   '', data_fn_train_ds  }, ...
%         { labels_fn_train, '', lbl_fn_train_ds }, ...
%         { mask_fn_train,   '', mask_fn_train_ds }, ...
%         data_fn_train );
%     ds.affinity_edges = [];
%     
%     dm.ds = ds;
%     this = dm;
%     
%     [~,prefix,~] = fileparts( lbl_fn_train_ds );
%     fnout = sprintf('%s/%s_infer_%s.h5', basedir, prefix, datestr(now,30) );
%     
%     dm.infer([],[],...
%         fnout, ...
%         [],[],[],[20 20 20],[], 3);
%     
% end

%%

if( ~ exist( data_fn_train_ds, 'file' ) )
    disp('writing data');
    im = read_image_stack( dm.ds.data_fn{1} );
    msk = read_image_stack( dm.ds.mask_fn{1} );
    lbl = read_image_stack( dm.ds.labels_fn{1} );
    
    im_ds  = reduceEffectiveResolution( im, dsFactor );
    msk_ds = reduceEffectiveResolution( msk, dsFactor );
    lbl_ds = reduceEffectiveResolution( lbl, dsFactor );
    
    writeH5( data_fn_train_ds, im_ds, 'single', [10 10 10]);
    writeH5( mask_fn_train_ds, msk_ds, 'uint8', [10 10 10]);
    writeH5( lbl_fn_train_ds,  lbl_ds, 'uint8', [10 10 10]);
end

%%

cube_size = 1;

x_test = 53;
y_test = 65;
z_test = 81;

xs = [ x_test x_test ];
ys = [ y_test y_test ];
zs = [ z_test z_test ];

data_size = dm.ds.get_data_size(flag,1);

flag = 3;

% dm.dds.dcs.dps = dm.dds.dcs.dps(2);
% dm.svm_normalization = 0;

feats = dawmr_feat_comp( dm, xs, ys, zs, data_size, flag, cube_size );

size( feats ) 
     

%% open the true and new results

% savedRes = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0008_crackSegDawmr_test/exp0008_crackSegDawmr_test_sub1_infer.h5';
% 
% seg = read_image_stack( savedRes );
% segnew = read_image_stack( fnout );
% 
% %%
% 
% % BrowseComponents('iii', seg, segnew, (seg - segnew))

%%

% savedRes = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0008_crackSegDawmr_test/exp0008_crackSegDawmr_test_sub1_infer.h5';
% seg = read_image_stack( savedRes );
% 
% dsList = [ 3 7 9 ];
% for i = 1:length(dsList)
%      
%     match = sprintf('%s/closeup_lbls_ds4_lr%d_infer_*.h5', basedir, dsList(i) )
%     a = dir(match);
%     fn = fullfile( basedir, a.name );
%     segds{i} = read_image_stack( fn );
%     
% end

%%

BrowseComponents('iiii', seg, segds{1}, segds{2}, segds{3} );