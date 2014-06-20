function out = convZ( im, ker, start )
% Usage:
%   out = convZ( im, ker, start )

sz = size( im );
len = length(ker);

imWin = zeros( [sz(1:2), len] );
rng = start : start+len - 1;
ii = (rng > 0 & rng <= sz(3));
rng = rng( ii );


imWin(:,:,ii) = im( :, :, rng );
kerRep = repmat( permute( ker, [ 1 3 2] ), sz(1), sz(2));

% squeeze(imWin)

out = sum( (imWin .* kerRep), 3);