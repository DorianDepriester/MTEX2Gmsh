mtexdata titanium
grains=calcGrains(ebsd);
grains=cond_smooth(grains);
G=gmshGeo(grains);
savegeo(G,'titanium.geo','ElementSize',10,'gradient',5,'ElementType','Tri');
