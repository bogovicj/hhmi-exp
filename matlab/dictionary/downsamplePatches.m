function [ patchesDs ] = downsamplePatches( patches, patchSize, dsFactor )
% Usage:
%   [ patches_ds ] = downsamplePatches( patches, dsFactor )
%
% Patches are stored in columns of the input patches array

[~, numPatches] = size( patches );

for i = 1:numPatches
    
    thisPatch = reshape( patches(:,i), patchSize);
    thisPatchDs = downSampleGaussianImglib( thisPatch, dsFactor );
    
    if( i == 1 )
       patchesDs = zeros( numel(thisPatchDs), numPatches ); 
    end
    
    patchesDs(:,i) = thisPatchDs(:);
    
end
