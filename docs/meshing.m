%%
mtexdata small
ebsd = ebsd('indexed');
grains = calcGrains(ebsd);
G=gmshGeo(grains);

%% Basic use
% The geometry described by the object |G| can be meshed using the Gmsh
% software as follows:
mesh(G,'default.msh')

%%
% The above command results in a 1 element thick mesh, consisting in
% linear wedge elements (6-node 3D elements. The element size is (roughly) 
% equal to the EBSD resolution.

%% Constant element size
% The default element size can be set as follows:
mesh(G,'constant_elmtSize.msh','ElementSize',50)
%%
% The resulting mesh cannot be (easily) displayed in MATLAB. Thus, the
% following illustrates the geometry when opening the mesh file with Gmsh:
% 
% <<msh_cst.png>>
% 

%%
% The unit here is the same as the EBSD map (ususally µm).

%% Size gradient
% Let $s(\mathbf{x})$ be the local element size at coordinates $\mathbf{x}$.
% The element size can be set as an increasing distance from the grains
% boundary such that:
%
% $$s(\mathbf{x})=s_0+kd_{GB}(\mathbf{x})$$
%
% with:
%
% * $s_0$ the element size at grain boundaries
% * $d_{GB}(\mathbf{x})$ the euclidean distance from the closest grain
% boundary.
%
% This can be done with the following command:
% 
%   mesh(G,meshFile,'ElementSize',s0,'gradient',k);
% 
%%
% 
% <<msh_gradient.png>>
% 


%% Element size Depending on the curvature of grain boundaries
% The local curvature of grain boundaries can be used to set the element
% size. For instance, the following command use 5 nodes to describe a full 
% circle:
mesh(G,'curvature.msh','Curvature',5);
%%
% 
% <<msh_curvature.png>>
% 


%% Element type
% The default element type for meshing is linear wedge. It can be
% changed to brick element
mesh(G,'brick.msh','ElementType','Brick');
%%
% 
% <<msh_brick.png>>
% 


%%
% or tetrahedrons
mesh(G,'tet.msh','ElementType','Tet');

%%
% If you wants to work in 2D only, use triangular (|Tri|) or quandrangular
% (|Quad|) elements instead:
mesh(G,'Tri.msh','ElementType','Tri');
mesh(G,'Quad.msh','ElementType','Quad');

%% Dump the geometry in an ASCII file
% The geometry can be exported into a Gmsh-readable (and somehow 
% human-readable) format using the following command:
savegeo(G,'geometry.geo')

%% 
% <html><hr></html>
%
% <index.html Go back to documentation index>
