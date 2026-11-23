mtexdata twins
grains=calcGrains(ebsd).smoothBoundary;
G=gmshGeo(grains);
grains
mesh(G,'twins.msh','ElementSize',0.2,'gradient',0.5,'elementType','Hex');