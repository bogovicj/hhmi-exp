function [gif_fn] = zflythruGif( im, fn, addcmds, scale, gifdelay)
% zflythruGif produces an animated gif that sequentially displays each 
% z slice in the input 3d volume
%
% Usage:
%   [gif_fn] = zflythruGif( im, fn );
%
% Additional plot options can be specified with the optional input
% 'addcmds'
%   [gif_fn] = zflythruGif( im, fn, 'axis equal; axis off;');
%
% Inputs:
%   im       - the input 3d volume
%   fn       - path of the output filename (excluding the .gif extension)
%   addcmds  - string input of additional commands to be executed
%              after the plotting
%   scale    - [0,1] scale data before rendering image?
%   gifdelay - time per frame in 0.01 sec( default = 20 )
%
% Outputs:
%   gif_fn - the output gif file name

if(~exist('addcmds','var'))
    addcmds = 'axis equal; axis off; colormap gray;';
end

if( ~exist('gifdelay','var') || isempty(gifdelay) )
    gifdelay = 10;
end

if( ~exist('scale','var') || isempty(scale) )
    scale = 1;
end

% make the directory for each image
mkdir( fn );

sz = size( im );

for i = 1:sz(3)
    
    figure('color','w');
    if ( scale )
        imagesc( im(:,:,i) );
    else
        imagesc( im(:,:,i) );
    end
    
    if ( ~isempty( addcmds ) )
       eval( addcmds ); 
    end
    
    % write to file here
    export_fig(fn,sprintf('%s/im_%05d.jpg',fn,i),'-nocrop','-native');
    close all;

end

% temporarily modify LD_LIBRARY_PATH for ImageMagick
old_path = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH', '/usr/local/cuda/lib64:');

% make the gif
gif_fn = sprintf('%s.gif',fn);
cmd = sprintf('convert -delay %d -loop 0 %s/*.jpg %s', ...
    gifdelay, fn, gif_fn);
fprintf('%s\n', cmd);
system(cmd);