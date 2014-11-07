function h = imdisp3d( im, varargin )
% Usage:
%   h = imdisp3d( im )

imdisp( permute( im, [1 2 4 3]), 'border', 0.02, varargin{:} );
h = gcf;
