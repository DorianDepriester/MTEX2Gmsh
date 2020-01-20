function Install_MTEX2Gmsh

%% Check if MTEX is installed
if ~isappdata(0,'MTEXInstalled')
	error('MTEX is not installed. It can be downloaded <a href="https://mtex-toolbox.github.io/download">here</a>.')
end

%% Check if Gmsh is installed and reachable
if ispref('MTEX2Gmsh','gmsh_path')
	path_to_gmsh=getpref('MTEX2Gmsh','gmsh_path');
else
	path_to_gmsh='gmsh';	% Default path
end

[val,ver]=system(sprintf('"%s" -version',path_to_gmsh));	% Check if Gmsh is reachable
if val
	msg='Gmsh has not be found in the PATH of your system. Do you wish to locate the executable file now?';
	title='Gmsh not found';
	no='Not now';
else
	msg=sprintf('Gmsh v%s has been found in your system. Use it for meshing?',strtrim(ver));
	title='Gmsh found';
	no='Choose another version';
end
answer=questdlg(msg, title,'Yes',no,'Yes');
if val && ~strcmp(answer,no) || ~val && strcmp(answer,no)
	if isunix
		[file,folder]=uigetfile('','Locate the executable file for Gmsh');
	else
		[file,folder]=uigetfile('*.exe','Locate the executable file for Gmsh');
	end
	if ischar(file) && ischar(folder)
		path_to_gmsh=[folder file];	% If operation cancelled, keep the default value
	end
end
setpref('MTEX2Gmsh','gmsh_path',path_to_gmsh)

%% Add MTEX2Gmsh to the path
fp=fileparts(mfilename('fullpath'));	% Path to the present file
if isempty(strfind(path,fp))
	addpath(fp);
end

%% Allows to install it permanently
answer = questdlg('Installation complete. Would you like to install MTEX2Gmsh permanently?', 'Install MTEX2Gmsh', 'Yes','No thanks','No thanks');
if strcmp(answer,'Yes')
	savepath;
end