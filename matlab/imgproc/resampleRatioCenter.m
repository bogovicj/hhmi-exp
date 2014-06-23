function out = resampleRatioCenter( im, ratio )
% Usage:
%   out = resampleRatioCenter( im, ratio )

%% test data 
% ratio = [1 1 2];
% im = rand( 13, 13, 13);

%% setup 

szin  = size( im );
szout = ceil(szin ./ ratio);

% halfWidthIn = (szin - 1) ./ 2

ctr_in = ceil(szin./2);
% ctr_out = ceil(szout./2);

kerHalfWidth = 1 + (ratio-1)./2;

%% do the work 

out = zeros( szout );

kz = [];
kz = resampCenterConvKernel( ratio(3), kz );

zc = ctr_in(3);

off = 0;
while ( ctr_in(3) + off * ratio(3) - kerHalfWidth(3) < szin(3))

    if( off == 0 )
        dd = 1;
    else
        dd = [ -1 1 ];
    end
       
    for d = dd
        zc = ctr_in(3) + d * off * ratio(3);
        zout = (zc + 1)./ratio(3);
        
        bi = zc - (length(kz) - 1)/2;
        out(:,:,zout) = convZ( im, kz, bi);
    end
    
    off = off + 1;
   
end

