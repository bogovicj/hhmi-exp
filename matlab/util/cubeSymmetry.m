function T = cubeSymmetry( i )
% T = cubeSymmetry( i )
% where i \in i = 32*n + 16*m + l;
% n in [ 0  1 ]
% m in [ 0  1 ]
% l in [ 0 11 ]
% 
% See:
% http://math.stackexchange.com/questions/78573/what-is-a-natural-way-to-enumerate-the-symmetries-of-a-cube 
%
% John Bogovic
% August 2014

if( ~exist('i','var') || isempty( i ))
    Tlist = cell(60,1); 
    n = 1;
    for i = 0:1
        for j = 0:1
            for k =0:11
                v = 32*i + 16*j + k;
                Tlist{n} = cubeSymmetry( v );
                n = n + 1;
            end
        end
    end
    T = Tlist( ~cellfun( 'isempty', Tlist ));
    return;
end

% initialize with the identity
T = eye( 3 );

% are we a reflection?
if( i >= 32 )
    T = [0 1 0; 1 0 0; 0 0 1] * T; % swap x and y
end

% do we swap tetrahedrons?
if( mod( i, 32 ) > 15 )
    T = [0 1 0; 1 0 0; 0 0 -1] * T;
end

% in tetrahedral group, peel of 120-ness
switch( bitshift( mod( i, 16 ), -2 ) )
    case 0
        % nothing
    case 1
%         T = T([2 3 1], : );
        T = [0 1 0; 0 0 1; 1 0 0] * T;
    case 2
%         T = T([3 1 2], : );
        T = [0 0 1; 1 0 0; 0 1 0] * T;
    otherwise
        error( 'invalid');
end

% in Klein four group, peel of 180-ness
switch ( mod( i, 4 ) )
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

