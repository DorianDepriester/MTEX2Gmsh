clear all
close all

%% MTEX stuff
plotx2east
mtexdata small
ebsd = ebsd('indexed');									% Remove unindexed points
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd);	% Compute grains

% region=[33800 5266 300 300];
% condition = inpolygon(ebsd,region);
% %ebsd = ebsd(condition)

selected_grains = grains(grains.grainSize > 1);			% Outlier removal
ebsd = ebsd(selected_grains);
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd);	% Update grains

%% gmshGeo stuff
%grains=grains([4 6]);
G=gmshGeo(grains);	% Format data structure

figure
plot(G);			% Plot surfaces

%% Generate Gmsh-readable file and export the mesh
savegeo(G,'small.geo','ElementSize',100,'gradient',1);

% Export grain properties
exportGrainProps(G,'small.csv');
Gmsh('small.geo','inp')


