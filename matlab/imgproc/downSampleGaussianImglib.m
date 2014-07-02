function im_ds = downSampleGaussianImglib( im, factor, sourceSigma, targetSigma )
% Gaussian downsampling of an image. 
% See http://fiji.sc/wiki/index.php/Downsample
%
% Usage:
%   im_ds = downSampleGaussian( im, newSize, sourceSigma, targetSigma )
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
    sourceSigma = [ 0 0 0.5 ];
end

if( ~exist('targetSigma','var') || isempty(targetSigma) )
    targetSigma = [ 0 0 0.5 ];
end

% sz = size( im );
% ndims = length( sz );
% newSz =  round( sz / factor );

isImglib = isa( im, 'net.imglib2.img.Img' );

if( ~ isImglib )
    imil = toImglib( im );
    im_ds_il =  net.imglib2.util.Resampling.resampleGaussian( ...
            imil, imil.factory(), ...
            factor, sourceSigma, targetSigma );
    im_ds = net.imglib2.util.ImgOps.toFloatArray3d( im_ds_il );
else
    im_ds =  net.imglib2.util.Resampling.resampleGaussian( ...
            im, im.factory(), ...
            factor, sourceSigma, targetSigma );
end
