%% renderBoxes

%%

alpha = 1;

%%
[verts, faces] = unitBox();

xColor = [ 0.5 0.2 0.2 ];
yColor = [ 0.2 0.5 0.2 ];
zColor = [ 0.2 0.2 0.5 ];

v = [verts ones(8,1)];
Tbase = eye(4);

scales = [1 1 3];

[t1, t2] = ndgrid( 0:2, 0:2 );
tvec = [t1(:) t2(:)];

for i = 1:9 
    
    thisverts = verts;
    
    T = Tbase;
    T(1:3,1:3) = diag( scales);
    T(1:2,4) = tvec(i,:)
    
    thisVerts = T*v';
    
    patch( 'faces', faces, 'vertices', thisVerts(1:3,:)', ...
           'edgecolor', xColor, ...
           'facecolor', 'none',  'facealpha', alpha );

end
    
axis equal; axis off; view(3);