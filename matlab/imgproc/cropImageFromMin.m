function [cimg, paramout] = cropImageFromMin(img,thresh,pad)
% [cimg, paramout] = cropImageFromMin(img,thresh,pad)
%   thres - voxels with intensities strictly less than thresh are
%   cropped

sz = size(img);
[i,j,k] = ind2sub(sz,find(img>thresh));
cimg = [];

if(isempty(i))
    error('no valid pixels meet the threshold criterion');
end

% length(i)
% lenth(j)
% length(k)
mini = min(i)-pad;
if(mini<1)
    mini=1;
end

maxi = max(i)+pad;
if(maxi>sz(1))
    maxi=sz(1);
end

minj = min(j)-pad;
if(minj<1)
    minj=1;
end

maxj = max(j)+pad;
if(maxj>sz(2))
    maxj=sz(2);
end

mink = min(k)-pad;
if(mink<1)
    mink=1;
end

maxk = max(k)+pad;
if(maxk>sz(3))
    maxk=sz(3);
end
paramout = [mini; maxi; minj; maxj; mink; maxk];
cimg = img(paramout(1):paramout(2), paramout(3):paramout(4), paramout(5):paramout(6));
