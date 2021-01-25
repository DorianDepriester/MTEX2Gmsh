function savegeo(obj,outputFilePath,varargin)
%SAVEGEO dumps the geometry into an ASCII file, using the Gmsh syntax.
% SAVEGEO('file.geo') dumps the geometry in file named file.geo.
%
% SAVEGEO('file.geo','option1',val1,'option2',val2,...) set name-value
% arguments for meshing option. The available options are the same as for
% the MESH function.
%
%	See also mesh.
	[~,~,fext] = fileparts(outputFilePath);
	if strcmpi(fext,'.geo')
		filepath = outputFilePath;
	else
		filepath = [outputFilePath '.geo'];	%	Append the extension if missing or wrong
	end	
	mesh(obj, filepath,varargin{:});		% Run the mesh command, except it doesn't mesh...
end