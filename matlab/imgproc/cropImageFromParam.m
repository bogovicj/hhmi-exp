function [cimg, paramout] = cropImageFromParam(img,param)
% [cimg, paramout] = cropImageFromParam(img,param)
paramout = param;
cimg = img(param(1):param(2), param(3):param(4), param(5):param(6));
