mtexdata titanium
grains=calcGrains(ebsd);
grains=cond_smooth(grains);
G=gmshGeo(grains);
mesh(G,'titanium_medium.msh','elementSize',40,'medium',[2000 2000 200],'mediumElementSize',180);
