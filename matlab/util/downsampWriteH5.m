function success = downsampWriteH5( srcfn, destfn, factor, op )
% Usage
%   success = downsampWriteH5( srcfn, destfn, factor, op )

success = 0;

try
    
    tmp = read_image_stack( srcfn );
    cls = class( tmp );
    sz = size( tmp );
    
    nc = sz(4);
    for c = 1:nc
        vol = tmp(:,:,:,c);
        
        
        volil = toImglib( vol );
        clear vol
        volilds = downSampleGaussianImglib( volil, factor );
        clear volil
        volDsC = net.imglib2.util.ImgOps.toFloatArray3d( volilds );
        
        if ( c == 1 )
            volDs = zeros( [ size(volDsC) nc ], cls );  
        end
        
        volDs(:,:,:,c) = volDsC;
        
        clear volilds volDsC
        
    end
    writeH5( destfn, volDs, cls );
    clear volDs
        
catch e
    e
end

success = 1;


%% OLD CODE
%     if( strcmpi( op, 'or' ) )
%         
%         vol = tmp(:,:,:,1);
%         for i = 2:sz(4)
%             vol = vol | tmp(:,:,:,2);
%         end
%         
%         
%     elseif( strcmpi( op, 'avg' ) )
%         vol = mean( tmp, 4 );
%         
%         
%     else
%         vol = tmp;
%        
%     end