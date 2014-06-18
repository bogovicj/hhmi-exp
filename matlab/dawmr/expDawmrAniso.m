%% expDawmrAniso
% experiment with dawmr clustering / segmentation using aniso data


%% gen synthetic data

% basedir = '/groups/saalfeld/home/bogovicj/projects/aniso/toys/grad';
% 
% sz  = [64 64 64 ];  % image size
% res = [ 1 1 1 ];    % resolution
% 
% t = [ 0.2 0.8 0.0];
% t = t./ sqrt( norm(t));
% [x,y,z] = meshgrid( 1:sz(1), 1:sz(2), 1:sz(3) );
% 
% im = t(1).* x + t(2).*y;
% clear x y z;
% 
% msk = ones( sz, 'uint8');
% lbl = msk;
% 
% im_fn = fullfile( basedir, 'im.h5');
% msk_fn = fullfile( basedir, 'msk.h5');
% lbl_fn = fullfile( basedir, 'lbl.h5');
% if( ~ exist( im_fn, 'file' ) )
%     sprintf('writing files');
%     writeH5( im_fn,   im, 'single', [10 10 10]);
%     writeH5( msk_fn, msk, 'uint8' , [10 10 10]);
%     writeH5( lbl_fn, lbl, 'uint8' , [10 10 10]);
% end
%% open an older model

savedDatFn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0008_crackSegDawmr_test/exp0008_crackSegDawmr_test_results.mat';

load( savedDatFn );

%% learn features

ds = dawmr_set( { data_fn_train,   '', data_fn_train  }, ...
                { labels_fn_train, '', labels_fn_train }, ...
                { mask_fn_train,   '', mask_fn_train }, ...
                data_fn_train );
ds.affinity_edges = [];

dm.ds = ds;
this = dm;

%% open the results

savedRes = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0008_crackSegDawmr_test/exp0008_crackSegDawmr_test_sub1_infer.h5';
seg = read_image_stack( savedRes );

%% grab a subset of the data 

x_test = 53;

nCorrect = 0;
nWrong = 0;

for x_test = 25:87
    

y_test = 65;
z_test = 81;
basedir = sprintf('/groups/saalfeld/home/bogovicj/projects/aniso/trn_sub_%d-%d-%d', x_test, y_test, z_test);
mkdir( basedir );


%%

% im  = read_image_stack( ds.data_fn{1} );
% msk = read_image_stack( ds.mask_fn{1} );
% lbl = read_image_stack( ds.labels_fn{1} );

% 
% im_test = im(   x_test-20 : x_test+20, ...
%                 y_test-20 : y_test+20, ...
%                 z_test-20 : z_test+20 );
% 
% lb_test = lbl(  x_test-20 : x_test+20, ...
%                 y_test-20 : y_test+20, ...
%                 z_test-20 : z_test+20 );
% 
% mk_test = msk(  x_test-20 : x_test+20, ...
%                 y_test-20 : y_test+20, ...
%                 z_test-20 : z_test+20 );
% 
% im_fn  = fullfile( basedir, 'im.h5');
% msk_fn = fullfile( basedir, 'msk.h5');
% lbl_fn = fullfile( basedir, 'lbl.h5');
% 
% if( ~ exist( im_fn, 'file' ) )
%     sprintf('writing files');
%     writeH5( im_fn,  im_test, 'single', [10 10 10]);
%     writeH5( lbl_fn, lb_test, 'uint8' , [10 10 10]);
%     writeH5( msk_fn, mk_test, 'uint8' , [10 10 10]);
% end


res_test = seg ( x_test, y_test, z_test )


%% run inference

% fnout =  sprintf('%s/infer_%s.h5', basedir, date );
% dm.infer([],[],...
%          fnout, ...
%          [],[],[],[20 20 20],[], 3);
     
%% look at result of inference

% seg_img = read_image_stack( fnout );
% BrowseComponents('ii', im, seg_img)

%% test out feature comp

xs = [ x_test x_test ];
ys = [ y_test y_test ];
zs = [ z_test z_test ];

temp_dir = fullfile( basedir, 'tmp');
if( ~exist( temp_dir, 'dir' ) )
   mkdir( temp_dir ); 
end

flag = 3;
cube_size = 20;
use_sparse = 0;
n_recur    = 0;
objInfer_fn = sprintf('%s/infer_dat.mat',temp_dir);

dm_fn  = sprintf('%s/dm_%s.mat', temp_dir, ...
    datestr(now,30));
save_info = dm.end_classifier.make_barebones();
if( ~exist( dm_fn, 'file' ) )
    save(dm_fn, 'this', '-v7.3');
end

full_data_size = dm.ds.get_data_size(flag,1);

dawmr_infer(dm_fn, xs, ys, zs, full_data_size, objInfer_fn, ...
                  flag, cube_size, use_sparse, n_recur);
              
%% check results

res = load( objInfer_fn );

if( abs(res.vals_pd_mat - res_test) > 0.0001 )
   fprintf( 'ERROR: expected %f \t got %f\n',  res_test, res.vals_pd_mat);
   
    nWrong = nWrong + 1;
else
   fprintf( 'CORRECT!\n' );  
   nCorrect = nCorrect + 1;
end

end 


nCorrect 
nWrong