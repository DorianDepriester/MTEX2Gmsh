mtexdata twins
grains = calcGrains(ebsd);
grains_smooth = grains.smoothBoundary;
G = gmshGeo(grains_smooth);
mesh(G, 'twins.msh', 'ElementSize', 0.2, 'gradient', 0.5, 'elementType', 'Hex');