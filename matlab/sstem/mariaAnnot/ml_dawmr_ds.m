function ds = ml_dawmr_ds( training, test )
% ds = ml_dawmr_ds( training, test )
%
% generate dawmr_set from maria luisa's annotations
%
% test is allowed to be empty

datdir = '/groups/saalfeld/home/bogovicj/projects/davi/mariaAnnot';
id2pre = {
    '140718_elastic_registered_AL_Maria-Luisa_ground_truth_v04', % 1
    '140721_elastic_registered_AL_Maria-Luisa_ground_truth_v03'  % 2
    };

if( ischar( training))
    trn = str2double( training );
else
    trn = training;
end

if( ~exist('test','var'))
    test = [];
end

if( ischar( test))
    tst = str2double( test );
else
    tst = test;
end

train_pre = id2pre{ trn };

data_fn_train   = sprintf('%s/%s/%s_im.h5',  datdir, num2str(trn), train_pre );
labels_fn_train = sprintf('%s/%s/%s_lb.h5',  datdir, num2str(trn),train_pre );
mask_fn_train   = sprintf('%s/%s/%s_msk.h5', datdir, num2str(trn),train_pre );

data_fn_test = '';
labels_fn_test = '';
mask_fn_test = '';

if( ~isempty( test ))
    test_pre = id2pre{ tst };
    data_fn_test   = sprintf('%s/%s_im.h5',  datdir, num2str(tst), test_pre );
    labels_fn_test = sprintf('%s/%s_lb.h5',  datdir, num2str(tst), test_pre );
    mask_fn_test   = sprintf('%s/%s_msk.h5', datdir, num2str(tst), test_pre );
else
%     data_fn_test   = data_fn_train;
%     labels_fn_test = labels_fn_train;
%     mask_fn_test   = mask_fn_train;
end

ds = dawmr_set( { data_fn_train,   '', data_fn_test   }, ...
                { labels_fn_train, '', labels_fn_test }, ...
                { mask_fn_train,   '', mask_fn_test   }, ...
                data_fn_train );
            
ds.affinity_edges = [];

