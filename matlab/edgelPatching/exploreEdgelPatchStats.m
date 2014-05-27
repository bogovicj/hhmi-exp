% exploring matching features for crack edgel patches
%
% run_script('exploreEdgelPatchStats', 'Visualize 2d harris edges');
% run_script('exploreEdgelPatchStats', 'Visualize 2d harris edges - quality 0.0001');
%
% run_script('exploreEdgelPatchStats', 'Visualize 2d FAST corners - quality 0.0001');

%% setup
global SAVEPATH;

ds = 4;
imf = sprintf('/groups/jain/home/bogovicj/projects/crackSegmentation/groundTruth/closeup/img_ds%d.tif', ds)
msf = sprintf('/groups/jain/home/bogovicj/projects/crackSegmentation/groundTruth/closeup/labels_interp_smooth_ds%d.tif', ds)


i_info = imfinfo(imf);
m_info = imfinfo(msf);

N = numel(i_info);
if( N ~= numel( m_info ) )
    error('Mask and Image are different sizes');
    %    return;
end

tmp = imread(imf);
I = zeros( [ size(tmp) N ]);
M = zeros( [ size(tmp) N ]);
for i=1:N
    I(:,:,i) = imread( imf, i, 'Info', i_info );
    M(:,:,i) = imread( msf, i, 'Info', m_info );
end



%% try harris

% minQuality = 0.0001; % decrease to get more points
% 
% % del = 5;
% del = 1;
% for i = 1:del:N
%     
%     slc = I(:,:,i)./255;
%     msk = M(:,:,i)./255;
%     mskDist = bwdist( ( msk > 0.5 ) );
%     
%     pts = detectHarrisFeatures( slc, 'MinQuality', minQuality );
%     
%     ptDists = mskDist( sub2ind( size(mskDist), ...
%                         round(pts.Location(:,1)), ...
%                         round(pts.Location(:,2))) );
%     
%     imagesc(slc); hold on; colormap('gray'); axis equal; axis off;
%     plot( pts ( ptDists < 6 ) );
%     
%     fn = fullfile( SAVEPATH, sprintf('slice_%04d.png',i));
%     fn
%     pts
%     
%     export_fig(fn, '-m2', '-q101');
% 
% %     % for debug
% %     pause;
% %     break; 
%        
%     close all;
% end

%% try fast

% minQuality = 0.0001; % decrease to get more points
% 
% % del = 5;
% del = 1;
% for i = 1:del:N
%     
%     slc = I(:,:,i)./255;
%     msk = M(:,:,i)./255;
%     mskDist = bwdist( ( msk > 0.5 ) );
%     
%     pts = detectFASTFeatures( slc, 'MinQuality', minQuality );
%     
%     ptDists = mskDist( sub2ind( size(mskDist), ...
%                         round(pts.Location(:,1)), ...
%                         round(pts.Location(:,2))) );
%     
%     imagesc(slc); hold on; colormap('gray'); axis equal; axis off;
%     plot( pts ( ptDists < 6 ) );
%     
%     fn = fullfile( SAVEPATH, sprintf('slice_%04d.png',i));
%     export_fig(fn, '-m2', '-q101');
% 
% %     % for debug
% %     pause;
% %     break; 
%        
%     close all;
% end


%% try surf

% del = 5;
% for i = 1:5:N
%     slc = I(:,:,i)./255;
%     
%     surf_points = detectSURFFeatures(slc );
% %     surf_points = detectSURFFeatures(slc,'MetricThreshold', 200);
%     
%     imagesc(slc); hold on; colormap('gray'); axis equal;
%     plot(surf_points.selectStrongest(10));
%     %    break;
%     pause;
%     clf;
% end
% 
% 
% %% try mser
% 
% del = 5;
% for i = 1:5:N
%     slc = I(:,:,i)./255;
%     
%     mser_reg = detectMSERFeatures(slc);
%     
%     imagesc(slc); hold on; colormap('gray'); axis equal;
%     plot(mser_reg(1:10),'showPixelList', true);
%     break;
%     %     pause;
%     %     clf;
% end

