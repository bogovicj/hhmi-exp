% analyze features

%% setup
featureFn = '/groups/jain/home/bogovicj/exp/saved_exp/exp0657_extractEdgelFeatures/Freak_features.mat';

load( featureFn );

Nslc = numel(i_info);
if( Nslc ~= numel( m_info ) )
    error('Mask and Image are different sizes');
    %    return;
end

tmp = imread(imf);
I = zeros( [ size(tmp) Nslc ]);
M = zeros( [ size(tmp) Nslc ]);
for i=1:Nslc
    I(:,:,i) = imread( imf, i, 'Info', i_info );
    M(:,:,i) = imread( msf, i, 'Info', m_info );
end

clear tmp;

%%
testFeature = features{end};
numCands = length(features) - 1;
metrics = zeros( numCands, 1);

for i = 1:numCands
    [matchIdx, matchMetric] = matchFeatures( testFeature, features{i}, 'MatchThreshold', 100 );
    metrics(i) = matchMetric;
end

hist( metrics );

%% 

figure;
testPtRnd = round(edgelPts(end,:));

matchIdx = find( metrics == min(metrics));


%% pick a match

ff = zeros( (length(features) - 1), size(features{1}.Features,2), 'uint8');
for i = 1:(length(features) - 1)
     ff(i,:) = features{i}.Features; 
end

[matchIdx, matchMetric] = matchFeatures( testFeature, binaryFeatures(ff), 'MatchThreshold', 100 );
matchIdx

%% nonsense

% % ff = zeros( (length(features) - 1), size(features{1}.Features,2), 'uint8');
% ff = zeros( 2, size(features{1}.Features,2), 'uint8');
% 
% % for i = 1:(length(features) - 1)
% for i = 1:2
%    ff(i,:) = features{i}.Features; 
% end
% 
% [matchIdx, matchMetric] = matchFeatures( testFeature, binaryFeatures(ff), 'MatchThreshold', 100, 'MaxRatio', 1 )