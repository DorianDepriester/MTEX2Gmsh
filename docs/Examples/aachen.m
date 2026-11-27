mtexdata aachen
grains=calcGrains(ebsd);
G=gmshGeo(grains.simplifyBoundary);
G.mesh('aachen.msh','ElementSize',0.7,'gradient',0.5)