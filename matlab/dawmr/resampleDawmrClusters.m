function resampleDawmrClusters( dm, ratio )
% resampleDawmrCluters
%
% replaces the contents of dm's clusters with a subsampled
% version

nds = length(dm.dds);

if( length(ratio) == 1 )
    r = [ 1 1 ratio ];
else
    r = ratio;
end

for di = 1:nds
    
    dd = dm.dds(di); ; % dawmr_data
    
    ndc = length( dd.dcs );
    
    for ci = 1:ndc
        dc = dd.dcs(ci); % dawmr_clustering
        
        cents = dc.centroids;
        for i = 1 : size(cents,1)
            cents_rs = reshape( cents(i,:), dc.patch_dim );
            cents_rs_sm = resampleRatioCenter( cents_rs, r );
            
            if ( i == 1 )
                patch_size_new = size( cents_rs_sm );
                cents_sm = zeros( size( cents,1), numel(cents_rs_sm));
            end
            cents_sm( i, : ) = cents_rs_sm(:)';
            
        end
        dc.patch_dim = patch_size_new;
        dc.centroids = cents_sm;
        
    end
    
end

