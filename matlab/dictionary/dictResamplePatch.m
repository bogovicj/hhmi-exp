function out = dictResamplePatch( in, dictionary, szout, resRatio )
% Usage:
% 	out = dictResamplePatch( in, dictionary, szout, resRatio );
%
% Inputs
%	in          - the input patch
%	dictionary  - the dictionary ( num_elements x num_variables)
%   dsszoutz    - size of the output ( must agree with size of dictionary patches )
%	resRatio    - resample ratio

n_dict = size( dictionary, 1);
sz = size( in );

numel_dict = size(dictionary,2);
numel_szout= prod(szout);
assert( (numel_szout == numel_dict), sprintf('dictonary size (%d) and size out (%d) do not agree', ...
    numel_dict, numel_szout));

% build the subsampled - dictionary
for i = 1:n_dict

   drs    = reshape( dictionary( i, :), szout );
   drssub = downSampleGaussianImglib( drs, resRatio );
  
   if( i == 1 )
		dict_sub = zeros( n_dict, numel( drssub ), class(dictionary));
   end
    
   dict_sub(i,:) = reshape( drssub, 1, [] );

end

in_rs = in(:);

alpha = dict_sub * in_rs;
out_rs = alpha' * dictionary;
out = reshape( out_rs, szout );
