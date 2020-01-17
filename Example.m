%% MTEX stuff 
plotx2east
mtexdata titanium
ebsd = ebsd('indexed');									% Remove unindexed points
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd);	% Compute grains

selected_grains = grains(grains.grainSize > 1);			% Outlier removal
ebsd = ebsd(selected_grains);
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd);	% Update grains
plot(grains)

%% gmshGeo stuff 
G=gmshGeo(grains);	% Format data structure
figure
plot(G);			% Plot the geometry

%% Generate the mesh and save it 
mesh(G,'small.inp','thickness',50,'elementType','Brick');

%% Export grain properties 
exportGrainProps(G,'small.csv');

