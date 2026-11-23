mtexdata titanium
grains=calcGrains(ebsd);
G=gmshGeo(grains.smoothBoundary);
mesh(G,'titanium-medium.msh','ElementSize',20,'medium',[2000 2000 200]);
