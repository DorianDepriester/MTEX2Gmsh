function mesh(obj,outputFilePath,varargin)
%MESH Mesh the geometry and export the mesh into the
%requested file.
%
%	MESH(obj,outputFile) meshes the geometry and saves the results
%	at the specified location. If the latter is only a filename
%	(no path), the mesh is saved in the current working directory. 
%	The element size is roughly equal to the EBSD resolution.
%
%	MESH(...,'ElementSize',s) results in element sizes equal to s.
%
%	MESH(...,'Thickness',t) sets an extrusion thickness equal to T 
%	(equal to element size by default).
%
%	MESH(...,'gradient',k) results in elements with size equal to
%	s+k*d (d being the distance from the nearest boundary and s the
%	default element size).
%
%	MESH(...,'ElementType',type) sets the element type used	for 
%	meshing. It can be:
%		-for 3D geometry:
%			-'Wedge' (default) for wedge elements,
%			-'Hex' for hexahedon elements,
%			-'Tet' or 'Tetrahedron' for tetrahedron elements,
%			-'HexOnly' for hexahedron elements only (no tet),		
%		-for 2D geometry:
%			-'Tri' or 'Triangular' for triangular elements,
%			-'Quad' or 'Quadrangular' for quadrangular elements,
%			-'QuadOnly' for quandrangular elements only (no
%			triangle).
%	Note that 'Quad' and 'Hex' lead to quad-dominated and
%	hex-dominated meshes, respectively, leaving some
%	wedge/triangles in the mesh; hence the 'HexOnly' and 'QuadOnly
%	options. 
%
%	MESH(...,'ElementOrder',order) sets the element order. The
%	default value is 1 (i.e. linear elements).
%
%	MESH(...,'Curvature',np) sets the element sizes to be computed
%	depending on the local curvature (np nodes per 2 pi). np==0 
%	disables this option (default).
%
%	MESH(...,'grainPrefix',str) defines the name for the element
%	sets corresponding to grains (Physical Volumes in Gmsh). E.g 
%	MESH(...,'grainPrefix','grain_') will create volumes named
%	'grain_1', 'grain_2' etc.
%	If the argument is empty, no prefix is given and the physical
%	volumes are just numbered as the grains.
%
%	MESH(...,'grainPrefix',str) defines the name for the element
%	sets corresponding to grains (Physical Volumes in Gmsh). E.g 
%	SAVEGEO(...,'grainPrefix','grain_') will create volumes named
%	'grain_1', 'grain_2' etc.
%	If the argument is empty, no prefix is given and the physical
%	volumes are just numbered as the grains.		
%
%	MESH(...,'medium',S) embeds the ROI inside a cuboid of size 
%	S=[dx dy dz]. The element size in the medium is	increasing with
%	increasing distance from the ROI. The mesh in the 	medium is 
%	composed of tetrahedron elements.
%
%	MESH(...,'medium',S,'mediumElementSize',value) sets the element
%	size at the corners of the medium to the given value.
%
%	See also savegeo.
	if ispref('MTEX2Gmsh','gmsh_path')
		path_to_gmsh=getpref('MTEX2Gmsh','gmsh_path');
	else
		path_to_gmsh='gmsh';	% Default path
	end
	[val,~]=system(sprintf('"%s" -version',path_to_gmsh));	% Check if Gmsh is reachable
	if val
		screenSize = get(0,'ScreenSize');	% Used for checking the -nodisplay option if off https://stackoverflow.com/a/6771356/12056867
		if isequal(screenSize(3:4),[1 1])	% -nodisplay option -> just cancel operation 
			error('Gmsh has not be found in the PATH of your system.')
		end
		answer=questdlg('Gmsh has not be found in the PATH of your system. You need to locate its executable file first.', 'Gmsh not found', 'Ok','Cancel','Cancel');
		if strcmp(answer,'Cancel')
			return
		end
		if isunix
			[file,path_to_gmsh]=uigetfile('','Locate the executable file for Gmsh');
		else
			[file,path_to_gmsh]=uigetfile('*.exe','Locate the executable file for Gmsh');
		end
		path_to_gmsh=[path_to_gmsh file];
		answer=questdlg('Would you like to save this option permanently?', 'Save preference', 'Yes','No','Yes');
		if strcmp(answer,'Yes')
			setpref('MTEX2Gmsh','gmsh_path',path_to_gmsh)
		end
	end			
	tmp_file=obj.savegeo(tempname,varargin{:});	%	Save the geometry into a temp file
	str=sprintf('"%s" "%s" -o "%s" -v 4 -3',path_to_gmsh,tmp_file,outputFilePath);
	system(str);
	delete(tmp_file)	%	delete temp file
end