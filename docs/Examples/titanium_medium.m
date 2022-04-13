mtexdata titanium
grains=calcGrains(ebsd);
grains=cond_smooth(grains);
G=gmshGeo(grains);
mesh(G,'titanium-medium.msh','ElementSize',20,'medium',[2000 2000 200]);
