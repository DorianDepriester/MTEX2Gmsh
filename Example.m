%% MTEX stuff
plotx2east
mtexdata small
ebsd = ebsd('indexed');									% Remove unindexed points
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd);	% Compute grains

selected_grains = grains(grains.grainSize > 1);			% Outlier removal
ebsd = ebsd(selected_grains);
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd);	% Update grains

%% gmshGeo stuff
G=gmshGeo(grains);	% Format data structure

figure
plot(G);			% Plot surfaces

%% Generate Gmsh-readable file and export the mesh
savegeo(G,'small.geo','thickness',50,'elementType','Brick');

%% Export grain properties
exportGrainProps(G,'small.csv');
Gmsh('small.geo','inp')


