function ebsd = ebsd_from_QuadratureOutputs(fname, CS)
%EBSD_FROM_QUADRATUREOUTPUTS reads quadrature outputs files from
%PRISMS-Plasticity and converts them into MTEX's ebsd data.
% 
% EBSD_FROM_QUADRATUREOUTPUTS('QuadratureOutputs.csv', CS) reads data from 
% the file named 'QuadratureOutputs.csv' and converts it to EBSD data,
% assuming a crystal symmetry given by CS.
%
% WARNING: Note that the quadrature outputs files given by 
% PRISMS-Plasticity can be very large. Thus this function may be slow are 
% even crash MATLAB. In this case, consider splitting those file using 
% external tools (see e.g. Python with Pandas or awk).

%
%	see also CrystalSymmetry

%% Preprocess data
T=readtable(fname); % Import results as table and keep only the columns of interest
Tmat=T{:,[2 5 6 8 9 10]};

% Convert Rodrigues components to Euler angles
rot=rotation.byRodrigues(vector3d(Tmat(:,4),Tmat(:,5),Tmat(:,6)));
eul=Euler(rot)/degree;	% use degrees!
Tmat(:,4:6)=eul;

% Remove NaN values
Tmat(any(isnan(Tmat),2),:)=[]; % In some cases (which ones?) some angles are NaN

%% Read data with MTEX
loader= loadHelper(Tmat,'ColumnNames', {'Phase' 'x' 'y' 'phi1','Phi','phi2'});
phase = loader.getColumnData('Phase');
q     = loader.getRotations();
opt = loader.getOptions('ignoreColumns','Phase');
ebsd = EBSD(q,phase,CS,opt);


end

