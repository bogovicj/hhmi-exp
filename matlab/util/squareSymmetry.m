function T = squareSymmetry( i )
% T = squareSymmetry( i )
% where i \in = [1:8];
% 
% John Bogovic
% September 2014

if( ~exist('i','var') || isempty( i ))
    Tlist = cell(8,1); 
    for i = 1:8
        Tlist{i} = squareSymmetry( i );
    end
    T = Tlist;
    return;
end

switch ( i )
    case 1
        T = eye( 2 );       % identity
    case 2
        T = [0 -1; 1 0];    % 90 deg rot
    case 3
        T = [-1 0; 0 -1];   % 180 deg rot
    case 4
        T = [0 1; -1 0];    % 270 deg rot
    case 5
        T = [-1 0; 0 1];    % flip x
    case 6
        T = [1 0; 0 -1];    % flip y
    case 7
        T = [0 1; 1 0];     % flip y=x diag 
    case 8
        T = [0 -1; -1 0];   % flip y=-x diag
    otherwise
        error('invalid')
end

