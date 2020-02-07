%% MTEX stuff
% Load EBSD data and reconstruct the grains
plotx2east
mtexdata titanium
ebsd = ebsd('indexed');									% Remove unindexed points
grains = calcGrains(ebsd);	% Compute grains
plot(grains)

%% gmshGeo stuff
% This converts grain2d data into Gmsh-like data
G=gmshGeo(grains);		% Format data structure
hold on
plot(G);				% Plot the geometry ontop of the grains
mesh(G,'small.inp');	% Export into Abaqus INP file

exportGrainProps(G,'small.csv');	% Export grain-related data into a CSV file