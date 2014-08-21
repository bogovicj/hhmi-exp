function T = cubeSymmetry( i )
% T = cubeSymmetry( i )
% where i \in {0:47}
% 
% See:
% http://math.stackexchange.com/questions/78573/what-is-a-natural-way-to-enumerate-the-symmetries-of-a-cube 

% initialize with the identity
T = eye( 3 );

% are we a reflection?
if bitand( i, 32 )
    T = T([2 1 3], : );
end

% do we swap tetrahedrons?
if bitand( i, 16 )
    T = T([2 1 3], : );
    T(3,3) = -1;
end

% in tetrahedral group, peel of 120-ness
switch( bitshift( bitand(i, 12), -2) )
    case 0
        % nothing
    case 1
        T = T([2 3 1], : );
    case 2
        T = T([3 1 2], : );
    otherwise
        error( 'invalid');
end

% in Klein four group, peel of 180-ness
switch ( bitand(i, 3) )
    case 0
        % nothing
    case 1
        T = [-1 0 0; 0 -1 0; 0 0 1] * T; % x y flip
    case 2
        T = [1 0 0; 0 -1 0; 0 0 -1] * T; % y z flip
    case 3
        T = [-1 0 0; 0 1 0; 0 0 -1] * T; % x z flip
    otherwise
        error('invalid')
end

%% test


% for i = 0:1
% for j = 0:1
% for k =0:11
% v = 32*i + 16*j + k;
% s{v+1} = cubeSymmetry(k);
% end
% end
% end
