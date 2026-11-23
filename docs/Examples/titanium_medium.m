mtexdata titanium
grains = calcGrains(ebsd);
grains_smooth = grains.smoothBoundary
G = gmshGeo(grains_smooth);
mesh(G, 'titanium-medium.msh', 'ElementSize', 20, 'medium', [2000 2000 200]);
