mtexdata twins
ebsd=ebsd('indexed')
grains=calcGrains(ebsd)
G=gmshGeo(cond_smooth(grains));
G=simplify(G);
mesh(G,'twins.msh','ElementSize',0.2,'gradient',0.5,'elementType','Brick');