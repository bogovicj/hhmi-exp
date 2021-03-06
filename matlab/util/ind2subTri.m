function [i,j] = ind2subTri(N,k)
% ind2subTri gives subscript indices (i,j) into a NxN trigular matrix stored as a
% vector, for example, generated by PDIST
%
% [i,j] = ind2subTri(N,k)
% where N is optional, as it can be calculated from k if it is of the
% proper form.  Might be slightly faster to provide N.

% if(length(D)>1)
%     M = length(D);
% else
%     N = D;
% end

% if(~exist('N','var') || isempty(N))
%     p = [1 -1 -2*N];
%     r = roots(p);
%     N = r(r>0);
%     if(mod(r,1)~=0)
%         error('Invalid vector - must have length N(N-1)/2, where N is an integer');
%     end
% end

%
% i d          
% 1 (N-1)
% 2 (N-1)+(N-2)     
% 3 (N-1)+(N-2)+(N-3)

i = ones(1,length(k));
j = ones(1,length(k));

for m = 1:length(k)
    
    d = (N-i(m));
    while(k(m) > d && i(m) < N)
        
        i(m) = i(m) + 1;
        d = d + (N - i(m));
    end
    
    % d - (N-i) gives the last index for the previous N-i long 'block'
    j(m) = i(m) + (k(m) - ( d - ( N-i(m) ) ));
    
end

%% an alternative method
% [i,j] = meshgrid( 1:N, 1:N);
% k = reshape( (i > j)', [], 1);
% 
% i = i(k);
% j = j(k);
% 
% i = i(n);
% j = j(n);


end