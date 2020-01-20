function Uninstall_MTEX2Gmsh
if ispref('MTEX2Gmsh')
	rmpref('MTEX2Gmsh');
end
fp=fileparts(mfilename('fullpath'));	% Path to the present file
rmpath(fp);
savepath;