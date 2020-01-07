function varargout=Gmsh(inputFilePath,outputFilePath)
%GMSH Execute Gmsh whithin Matlab.
%   GMSH(INPUT,OUTPUT) runs the Gmsh sofware using the geometry specified
%   in INPUT and exports the mesh in the specified OUTPUT file. Gmsh must 
%   belong to the PATH of your operating system.
%
%	The OUTPUT argument can be:
%		- file extension (eg. '.inp'). In this case, the output file
%		is at the same location as the input file and is named identically
%		with specified extension.
%		- file format (eg 'INP'). This is just a shortcut for the
%		corresponding file extension (eg 'INP' -> '.inp').
%		- full path (eg 'example.inp')
%
%   By default, the returned message from Gmsh is printed in the Matlab's 
%   command window.
%
%   s=GMSH(...) saves the message as a string.
%
%	See also gmshGeo/mesh gmshGeo/savegeo
	
	%% Remove file extension from input file
	[inputPath,inputFileName,~] = fileparts(inputFilePath);
	if ~isempty(inputPath)
		inputPath=[inputPath filesep];
	end
	
	%% Check if the file exists
	if ~exist(inputFilePath,'file')
		error([inputFilePath ' not found']);
		return
	end	
	
	%% Process the output file path
	[outputPath,outputFileName,ext]=fileparts(outputFilePath);
	if isempty(outputPath)	% No path specified: use the same location as input file
		if isempty(ext)		% Format is given as a simple string (eg 'INP')
			if isempty(outputFileName)
				error('The output argument can not be empty.')
			end
			outputFilePath=[inputPath inputFileName '.' lower(outputFileName)];
		else				% File extension is given per se
			if isempty(outputFileName)
				outputFileName=inputFileName;
			end
			outputFilePath=[inputPath outputFileName ext];
		end
	else
		if isempty(ext)
			error('File extension is needed when giving a full path.');
		end
	end
	
	%% Format and run the command line
	str=sprintf('gmsh "%s" -o "%s" -3',inputFilePath,outputFilePath);	
	if nargout==0
		system(str);												% Run Gmsh
	else
		[~,varargout{1}]=system(str);								% Run Gmsh and get the output message
	end
end

