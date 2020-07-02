setMTEXpref('generatingHelpMode',true); % Avoid some artefact (fix issue #5)
%% MTEX stuff
% In this example, EBSD are loaded from the dataset named 'small', provided
% along with the MTEX toolbox. Once loaded, they can used for
% reconstructing grains (of class |grain2d| ).
% 
mtexdata small
ebsd = ebsd('indexed');		% Remove unindexed points
grains = calcGrains(ebsd);	% Compute grains
plot(grains)

%% Compute the meshable geometry and plot it
% The command below converts the grain geometries into a Gmsh-like
% description of the domains:

G=gmshGeo(grains);

%%
% The variable G is an object of class |gmshGeo| . It fully describes the
% geometries of all grains and their intrisic properties.

%% Plot the geometry
% One can easily plot the geometry using the usual plot function. 
plot(G);

legend show

%% Mesh the geometry

mesh(G,'small.inp');

%% 
% This command runs Gmsh within Matlab and exports the mesh into the
% specified format (Abaqus INP here). Extensive tweaks are available for
% meshing, as detailed in <mesh.html the corresponding section> .

%% Export the grain properties
% In order to retrieve the grain properties when importing the mesh into
% the FEM code (for assigning materials and orientations), one can export
% the properties of each grain as Comma-separated values (CSV).
exportGrainProps(G,'small.csv');

%% 
% <html><hr></html>
%
% <index.html Go back to documentation index>