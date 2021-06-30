function exportGrainProps(obj, filename, varargin)
%EXPORTGRAINPROPS Exports grain properties as ASCII data in a CSV 
% file.
%
%	EXPORTGRAINPROPS(G,'filename') exports grain properties stored
%	in G in the ASCII file named 'filename'. By default, the dumped 
%	properties are (in that order):
%		- GrainID,
%		- Phase name,
%		- Euler angle phi1,
%		- Euler angle Phi,
%		- Euler angle phi2.
%
%	EXPORTGRAINPROPS(G,'filename','prop1','prop2',...) exports the
%	properties named 'prop1', 'prop2' and so on, in that order. The
%	available properties are:
%		- 'GrainID',
%		- 'Phase',
%		- 'PhaseID',
%		- 'Euler',
%		- 'Rodrigues'.
%
%	The 'phaseID' option will number all the unique phases and use
%	those indices instead of full phase names (as in 'Phase').
%
%	The 'Euler' option stands for all Euler angles 'phi1', 'Phi',
%	'phi2'.
%
%	The 'Rodrigues' option converts the Euler angles into the
%	Rodrigues coordinates (r_x, r_y, r_z).
%
%	The 'PRISMS' option makes the output file suitable for 
%	PRISMS-Plasticity [1]. 
%	It is an alias for ' 'GrainID','Rodrigues', 'PhaseID' '.
%
%	[1] M. Yaghoobi, S. Ganesan, S. Sundar, A. Lakshmanan, S. 
%	Rudraraju, J.E. Allison, V. Sundararaghavan, 
%	"PRISMS-Plasticity: An open-source crystal plasticity finite 
%	element software” Computational Materials Science 169 (2019). 
%
%	See also savegeo.
	data=obj.Grains;
	if isempty(varargin)
		flags={'GrainID','Phase','Euler'};
	elseif length(varargin)==1	&& strcmpi(varargin{1},'PRISMS')
		flags={'GrainID','Rodrigues','PhaseID'};
	else
		flags=varargin;
		valid_args=horzcat(data.Properties.VariableNames,{'PhaseID','Euler','Rodrigues'});
		if ismember('PRISMS',flags)
			error('''PRISMS'' cannot be used together with other arguments.')
		end				
		if ~all(ismember(flags,valid_args))
			error('Wrong property name. It can be any combination of:\n - ''%s''\n or ''PRISMS'' alone.', strjoin(valid_args,''',\n - '''));
		end
	end
	[a,b]=ismember({'Euler'},flags);
	if a	% User asks for Euler angles
		flags=horzcat(flags{1:(b-1)},{'phi1','Phi','phi2'},flags{(b+1):end});
	end
	[a,b]=ismember({'Rodrigues'},flags);
	if a	% User asks for Rodrigues rotation
		Rodrigues_var_names = {'r_x','r_y','r_z'};
		flags=horzcat(flags{1:(b-1)},Rodrigues_var_names,flags{(b+1):end});
		rot=rotation.byEuler(obj.Grains.phi1,obj.Grains.Phi,obj.Grains.phi2);
		rodr=-Rodrigues(rot); % see https://groups.google.com/g/prisms-cpfe-users/c/4B4gRXXa5U0/m/ZJyGmYJaAgAJ
		newcols = table(rodr.x,rodr.y,rodr.z,'VariableNames',Rodrigues_var_names);
		data = horzcat(data,newcols);
	end
	[a,b]=ismember({'PhaseID'},flags);
	if a	% User asks for phase ID
		flags=horzcat(flags{1:(b-1)},{'PhaseID'},flags{(b+1):end});
		phaseList=unique(data.Phase);
		[~,PhaseID]=ismember(data.Phase,phaseList);
		newcol = table(PhaseID);
		data = horzcat(data,newcol);
		fprintf('In %s, the phases are numbered as follows:\n',filename);
		for i=1:length(phaseList)
			fprintf('  %i: %s\n',i,phaseList{i});
		end
	end	
	writetable(data(:,flags),filename,'delimiter','\t')
end