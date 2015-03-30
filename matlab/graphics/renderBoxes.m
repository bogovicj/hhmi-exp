%% renderBoxes

%%

alpha = 0.7;
[verts, faces] = unitBox();

%%


xColor = [ 0.5 0.2 0.2 ];
yColor = [ 0.2 0.5 0.2 ];
zColor = [ 0.2 0.2 0.5 ];

v = [verts ones(8,1)];
Tbase = eye(4);

scales = [1 1 2];

[t1, t2] = ndgrid( 0:4, 0:4 );
tvec = [t1(:) t2(:)];
zoff = 3;

for i = 1:size(tvec,1)
    
    thisverts = verts;
    
    T = Tbase;
    T(1:3,1:3) = diag( scales);
    T(1:2,4) = tvec(i,:);
    T(3,4) = zoff;
    
    thisVerts = T*v';
    
    patch( 'faces', faces, 'vertices', thisVerts(1:3,:)', ...
           'edgecolor', 'k', ...
           'facecolor', xColor,  'facealpha', alpha );

end
    
axis equal; axis off; view(3);

%
rng = 0.5:8.5;
[px,py,pz] = ndgrid(rng,rng,rng);
hrpts = [ px(:), py(:), pz(:)];
ptargs = { 'facecolor', 'k', 'edgecolor', 'none' };
[sf,sv]=renderPointSpheres( hrpts, gcf, 0.1, [], ptargs );

% do the bounding box of the HR patch

Tbbox = Tbase;
Tbbox(1:3,1:3) = diag( [9 9 9] );
bboxv = Tbbox * v';

hold on;
patch( 'faces', faces, 'vertices', bboxv(1:3,:)', ...
       'edgecolor', 'k', 'facecolor', 'none' );

