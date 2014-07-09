function upInd = upsampleIndicator( szLo, factor )
% 
% 
% Usage:
%   upInd = upsampleIndicator( szLo, factor )


% szLo   = [ 5 5 3 ];
% factor = [ 1 1 2 ];

szHi = (szLo - 1) .* factor + 1;
[x,y,z] = meshgrid(1:szLo(1), 1:szLo(2), 1:szLo(3));

xHi = (x-1).*factor(1) + 1;
yHi = (y-1).*factor(2) + 1;
zHi = (z-1).*factor(3) + 1;

upInd = false( szHi );
upInd( xHi, yHi, zHi ) = true;
