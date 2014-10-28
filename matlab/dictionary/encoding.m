function feat = encoding( X, D, encoding_type, do_rev_pol, param )
% feat = encoding( X, D, encoding_type, do_rev_pol, param )
%
% In:
%   X - matrix( M x N )
%   D - matrix( M x F )
%
% Out:
%   feat - matrix( F x M )
%
% M - length of data vector
% N - number of data samples
% F - size of dictionary
%
% encoding_type 
%   'sc'   - sparse coding
%   'st'   - soft thresholding
%   'best' - best

% alpha = D'*X; % F x N

if( ~exist('encoding_type','var') || isempty( encoding_type ))
    encoding_type = '';
end
if( ~exist('do_rev_pol','var') || isempty( do_rev_pol ))
    do_rev_pol = '';
end

switch encoding_type
    case 'sc'
        fprintf('sparse coding encoding...\n');
        feat = mexLasso(X, D, param);
    case 'st'
        fprintf('soft threshold encoding...\n');
        
        feat =  D'*X;
        fmn = mean(abs(feat),1);
        tmp =  bsxfun( @lt, feat, fmn) & bsxfun( @gt, feat, -fmn);
        feat(tmp) = 0;

    case 'best'
        feat =  D'*X;
        [~, fmi] = max( feat, [], 1);
        tmp = true( size( feat ));
        for i = 1:size(feat,2)
            tmp(fmi(i),i) = false;
        end
        feat(tmp) = 0;
        
    otherwise 
        disp('otherwise');
        
end

if( do_rev_pol )
    feat_raw = feat;
    feat = [ max(0, feat_raw) max(0, -feat_raw)];
end