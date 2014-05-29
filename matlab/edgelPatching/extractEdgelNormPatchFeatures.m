% exploring matching features for crack edgel patches
%
% dbstop if error; run_script('extractEdgelNormPatchFeatures', 'FREAK features ds 4');
% dbstop if error; run_script('extractEdgelNormPatchFeatures', 'FREAK features ds 4, testPt');
%
% dbstop if error; run_script('extractEdgelNormPatchFeatures', 'FREAK features ds 2, testPt');

%% setup
global SAVEPATH;
global DATDIR;

if ( strcmp( getenv('HOSTNAME'), 'jainlab-ws5.janelia.priv') )
   DATDIR = '/data-ssd1/john';
else
    DATDIR = '/groups/jain/home/bogovicj';
end

% patchDir = '/data-ssd1/john/projects/crackPatching/closeup/patches/matchExpsDepth_ds_4_31-31-19/test_8345';
% patchDir = '/data-ssd1/john/projects/crackPatching/closeup/patches/matchExpsDepth_ds_2_45-45-19/test_35193';

% fList = dir(fullfile(patchDir,'*Patch*tif'));    
% N = length(fList)

%% setup / preliminaries

downSampleFactor = 4;
searchRadius = 100;
searchCount  = 2000;

winAvgRad = 0;

testPt = [ 173, 205, 408 ]; 
testPt = testPt./downSampleFactor;

imgfn = sprintf('%s/projects/crackSegmentation/groundTruth/closeup/img_ds%d.tif', DATDIR, downSampleFactor)
maskfn = sprintf('%s/projects/crackSegmentation/groundTruth/closeup/labels_interp_smooth_ds%d.tif', DATDIR, downSampleFactor)
img =  net.imglib2.img.ImagePlusAdapter.convertFloat( ij.IJ.openImage( java.lang.String(imgfn)) );
mask = net.imglib2.img.ImagePlusAdapter.convertFloat( ij.IJ.openImage( java.lang.String(maskfn)) );
		
patchSize = [45, 45, 19 ];

cc = net.imglib2.algorithms.crack.CrackCorrection( img, mask, patchSize );
cc.computeEdgels();
edgels = cc.getEdgels();
nEdgels = edgels.size();

em = net.imglib2.algorithms.crack.EdgelMatching(img, mask, patchSize);
em.debugDir = sprintf('%s/projects/crackPatching/cropEdgelMatch', DATDIR);
em.debugSuffix = sprintf('ds%d',downSampleFactor);
		
cc.edgelMatcher = em;
cc.edgelMatcher.setEdgelSearchRadius( searchRadius / (downSampleFactor) );

eidx = cc.edgelIdxNearest(testPt);
testEdgel = cc.getEdgels().get(eidx);

%% edgel matching

cc.edgelMatcher.setEdgels( edgels );
matches = cc.edgelMatcher.candidateEdgels( testEdgel );
cc.edgelMatcher.filterEdgelsByNormal( testEdgel, matches );
nMatches = matches.size();
nMatches

% % for debug
% nMatches = 20;

%%
factory = net.imglib2.img.array.ArrayImgFactory();
tmp = factory.create( patchSize, img.firstElement());

midpt = (patchSize-1)/2 + 1;
midz  = midpt(3);

rl = org.apache.log4j.LogManager.getRootLogger();
rl.setLevel(org.apache.log4j.Level.INFO);

allFeatures = zeros( nMatches, 64, 'uint8' );
edgelPts    = zeros( nMatches, 3);
testFeature = [];
for i = -1:(nMatches-1)
% for i = -1
    
    if( mod (i,10) == 0 )
        fprintf('working on edgel %d of %d\n', i, nMatches);
    end

    if( i == -1 )
        e = testEdgel;
    else
        e = matches.get(i);
        edgelPts(i+1,1) = e.getDoublePosition(0);
        edgelPts(i+1,2) = e.getDoublePosition(1);
        edgelPts(i+1,3) = e.getDoublePosition(2);
    end
    
    v = net.imglib2.algorithms.edge.EdgelTools.edgelToView( e, img, patchSize );
    net.imglib2.util.ImgOps.copyInto( v, tmp );
    imar = net.imglib2.util.ImgOps.toDoubleArray3d( tmp );
    
    slc = mean(imar( :, :, midz-winAvgRad : midz + winAvgRad ), ...
                3);
    slc = slc./255;
    
    [feats, pts] = extractFeatures(slc, midpt(1:2), 'Method', 'FREAK');
    if ( i == -1 )
        testFeature  = feats.Features;
    else
        allFeatures((i+1), :) = feats.Features;
    end
    
    
end

%% analysis 

sz = [ img.dimension(0), img.dimension(1), img.dimension(2) ];
metrics = zeros( nMatches, 1);
metricImg = zeros(sz ); % tiff image visualizing output

for i = 1:nMatches
    
%     if ( isempty(features{i}) )
%         continue;
%     end
       
    [matchIdx, matchMetric] = matchFeatures( testFeature, allFeatures(i,:), 'MatchThreshold', 100 );
    metrics(i) = matchMetric;
    ptRnd = round( edgelPts(i,:));
    metricImg( ptRnd(2), ptRnd(1), ptRnd(3)) = metrics(i);
    
end

%% write tiff of metric info

metricfn = fullfile( SAVEPATH,'metricImg.tif');
net.imglib2.util.ImgOps.writeToTiff( metricImg, metricfn );


%% clean up
clear imar slc cc edgels matches testEdgel tmp v metricImg mask img e em feats factory rl

%%
% for debug
% SAVEPATH = '/data-ssd1/john/projects/crackPatching/closeup' 

fn = fullfile( SAVEPATH, 'results.mat' );
save( fn );

%% some older code

% midpt = [];
% features = cell(N,1);
% hasResults = false( N, 1 );
% for i=1:N
%     
%     f = fullfile( patchDir, fList(i).name);
%     img = imread( f );
%     
%     if( i == 1)
%         midpt = (size(img)-1)/2 + 1;
%     elseif ( mod(i,10) == 0 )
%         fprintf('patch %d of %d\n', i, N);
%     end
%     
%     [feats, pts] = extractFeatures(img, midpt, 'Method', 'FREAK');
%     feats
%     if( feats.NumFeatures > 0 )
%         hasResults( i ) = true;
%         features{i} = feats;
%     end
%     pause;
% end
