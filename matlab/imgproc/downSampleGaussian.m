function im_ds = downSampleGaussian( im, factor, sourceSigma, targetSigma )
% Gaussian downsampling of an image. 
% See http://fiji.sc/wiki/index.php/Downsample
%
% Usage:
%   im_ds = downSampleGaussian( im, factor, sourceSigma, targetSigma )
%
% newSize - a scalar specifying the downsampling factor.
%
% Motivation:
% Sound downsampling of an image requires the elimination of image frequencies
% higher than half the sampling frequency in the result image (see the
% Nyquist-Shannon sampling theorem).  The exclusive tool for this is Gaussian
% convolution.
%
% This script calculates the required Gaussian kernel for a given target size,
% smoothes the image and resamples it.
%
% Furthermore, you can define the "intrinsic" Gaussian kernel of the source and
% target images.  An optimal sampler is identified by sigma=0.5.  If your
% source image was blurred already, you may set a higher source sigma for a
% sharper result.  Setting target sigma to values smaller than 0.5 makes the
% result appear sharper and therefore eventually aliased.

if( ~exist('sourceSigma','var') || isempty(sourceSigma) )
    sourceSigma = 0.5;
end

if( ~exist('targetSigma','var') || isempty(targetSigma) )
    targetSigma = 0.5;
end

sz = size( im );
ndims = length( sz );
newSz =  round( sz / factor );

s = targetSigma * factor;
sig = sqrt( s * s - sourceSigma * sourceSigma );


if( ndims == 2 )
    im_ds = imresize( ...
        imfilter(im, fspecial('gaussian', [3 3], sig) ), ...
        newSz, 'box', 'antialiasing', 0 );
    % im_ds = imfilter(im, fspecial('gaussian', [3 3], sig) );
elseif( ndims == 3 )
   im_ds = smooth3( im, 'gaussian', sig ); 
end

         