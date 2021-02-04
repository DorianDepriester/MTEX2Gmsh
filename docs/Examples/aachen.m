mtexdata aachen
ebsd=ebsd('indexed');
grains=calcGrains(ebsd);
G=gmshGeo(grains);
G=simplify(G);
mesh(G,'aachen.msh','ElementSize',0.7,'gradient',0.5)