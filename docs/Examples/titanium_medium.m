mtexdata titanium
grains=calcGrains(ebsd);
grains=cond_smooth(grains);
G=gmshGeo(grains);
mesh(G,'titanium.geo','ElementSize',10,'gradient',5,'ElementType','Tri');
