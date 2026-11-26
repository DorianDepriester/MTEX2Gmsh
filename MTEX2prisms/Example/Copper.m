%% Example file for MTEX2prisms
% This script is an example illustrating how EBSD can be converted to a mesh
% suitable for PRISMS-Plasticity. It generates two files:
%	- the mesh
%	- the tabular file specifying the orientations of each grain, using the
%	Rodrigues convention.

%% Load EBSD data
mtexdata copper
% This is an example dataset provided by MTEX

%% Keep indexed point only
ebsd = ebsd('indexed');

%% Compute grains
[grains, ebsd] = calcGrains(ebsd);

%% Remove small grains
ebsd(grains(grainNumPixel(grains)<=20)) = [];
[grains, ebsd] = calcGrains(ebsd);

%% Smooth GBs
grains_smooth=grains.smoothBoundary;
plot(grains_smooth,grains_smooth.meanOrientation, 'noBoundary')

%% Move bottom-left corner to (0,0)
V=grainVertices(grains_smooth);
ori=min(V);
grains_smooth=setGrainVertices(grains_smooth,V-repmat(ori, size(V, 1), 1));

%% Convert data with MTEX2Gmsh
G=gmshGeo(grains_smooth);
hold on
plot(G)

%% Mesh geometry
G.mesh('Copper.msh','elementSize',5,'gradient',0.7,'ElementType','HexOnly','thickness',10,'grainPrefix','','verbosity', 0);
% Setting verbosity to 0 avoids thousands of warnings.
% You sometimes have to play around with the elements size options in order 
% to get sufficient mesh quality.

%% Export grain properties
G.exportGrainProps('orientations.txt','PRISMS')