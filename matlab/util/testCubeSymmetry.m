%% testCubeSymmetry

[x,y,z] = meshgrid( -1:1, -1:1, -1:1 );
pts = [ x(:), y(:), z(:) ];
size(pts)


clear s
s={};
for i = 0:1
    for j = 0:1
        for k =0:11
            v = 32*i + 16*j + k;
            s{v+1} = cubeSymmetry( v );
            fprintf('\n\n');
        end
    end
end

s = s( ~cellfun( 'isempty', s ));

ptXfm = [];
trs = [];

length( s )

k = 1;
for i = 1:length(s)
    trs(k,:) = s{i}(:);
    
    res = s{i} * pts';
    res = res';
    ptXfm(k,:) = res(:);
    k = k + 1;
    
end

%%

dt = pdist( trs );
nnz( dt == 0 )

dts = squareform( dt );
find( dts(1,:)==0 )

%%
D = pdist( ptXfm );
size(D)

nnz( D == 0 )