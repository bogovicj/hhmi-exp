function [ offset, funVal, xc, noverlap ] = matchByMaximalPMCC( template, image, maxCurv, requiredNumberOfOverlapPixels )
% [ offset, funVal ] = matchByMaximalPMCC( template, image, maxCurv ) 
%   template and image must be 2D 

if( ~exist('maxCurv','var') || isempty( maxCurv ))
    maxCurv = 100;
end

if( ~exist('requiredNumberOfOverlapPixels','var') || isempty( requiredNumberOfOverlapPixels ))
    requiredNumberOfOverlapPixels = 0.1 .* numel( template );
elseif ( requiredNumberOfOverlapPixels > 0  && ...
        requiredNumberOfOverlapPixels < 1 )
    
    requiredNumberOfOverlapPixels = requiredNumberOfOverlapPixels .* numel( template );
end

[ xc, noverlap] = normxcorr2_general( template, image, requiredNumberOfOverlapPixels );

midpt = (size(xc) - 1)./2;

out = javaArray( 'java.lang.Double', 2 );
funVal = mpicbg.ij.blockmatching.BlockMatching.quadraticFit( ...
    xc, midpt(1), midpt(2), maxCurv, out);

offset = zeros(2,1);
offset(1) = out(1);
offset(2) = out(2);

%% an old test
% % clc;
% % clear java; startup;
% 
% % a = zeros( 9 ); a( 1:4, : ) = 1;
% % b = zeros( 9 ); b( 1:5, : ) = 1;
% % a = a + 0.02 .* rand(9);
% % b = b + 0.02 .* rand(9);
% 
% a = zeros( 9 ); 
% % a = 0.2 .* rand( 9 );
% 
% a( 4, : ) = a( 4, : ) + 0.8 ;
% a( 5, : ) = a( 5, : ) + 1.2;
% a( 6, : ) = a( 6, : ) + 1.1;
% a( :, 4 ) = a( :, 4 ) - 0.2; 
% a( :, 6 ) = a( :, 6 ) - 0.2;
% a(4:6,4:6)
% 
% out = javaArray( 'java.lang.Double', 2 );
% 
% res = mpicbg.ij.blockmatching.BlockMatching.quadraticFit( ...
%     a, 4, 4, 100, out);
% 
% out
% res
