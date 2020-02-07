%% Prerequities
% In order to install MTEX2Gmsh in MATLAB, first check that: 
%
% * the <https://mtex-toolbox.github.io/ MTEX toolbox> for MATLAB is
% already installed,
% * <http://gmsh.info/ Gmsh> is installed on your computer.
% 

%% Install MTEX2Gmsh
% Just run the following command for installing the MTEX2Gmsh toolbox:
Install_MTEX2Gmsh


%%
% This function will check that the aforementioned depedencies are
% satisfied and eventually ask a few questions about the location of the 
% Gmsh executable file. This location will be stored persistently (if you 
% close MATLAB and re-open it, those informations will remain). Rerun the
% previous command if your configuration changes.

%% Uninstall MTEX2Gmsh
% If you have finished with this toolbox, you can uninstall it with:
Uninstall_MTEX2Gmsh
