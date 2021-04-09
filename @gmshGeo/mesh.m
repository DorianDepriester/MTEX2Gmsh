function fh=mesh(obj,filepath,varargin)
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
%	MESH(...,'periodic',axis) adds periodic condition on mesh along the
%	given axis. The axis can be 'X', 'Y' or 'both'. 'None' disables the
%	periodicity (default).
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
%	MESH(...,'medium',S) embeds the ROI inside a cuboid of size 
%	S=[dx dy dz]. The element size in the medium is	increasing with
%	increasing distance from the ROI. The mesh in the 	medium is 
%	composed of tetrahedron elements.
%
%	MESH(...,'LocalSize',S) sets local element size for specified grains. S
%	is a [N 2] matrix (with N the number of grains where local size is 
%	prescribed). The structure of S must be as follows:
%		S=[	grainID_1  size_in_grain1;
%			grainID_2  size_in_grain2;
%			...			...			]
%
%	MESH(...,'medium',S,'mediumElementSize',value) sets the element
%	size at the corners of the medium to the given value.
%
%	MESH(...,'verbosity',value) sets verbosity for Gmsh (0 for silent, 10
%	for max verbosity. Default: 4.
%
%	MESH(...,'partition',p) will create p partitions in the mesh for
%	parallel compouting.
%
%	See also savegeo, addSymmetry.

	%%	Parse optional parameters
	p = inputParser;
	addOptional(p,'ElementSize',0);
	addOptional(p,'thickness',0);
	addOptional(p,'gradient',0);
	addOptional(p,'ElementType','Wedge');
	addOptional(p,'ElementOrder',1);
	addOptional(p,'Curvature',0);
	addOptional(p,'Medium',[0 0 0]);
	addOptional(p,'MediumElementSize',0);
	addOptional(p,'grainPrefix','Grain_');
	addOptional(p,'verbosity',4);
	addOptional(p,'partition',0);
	addOptional(p,'periodic','none');
	addOptional(p,'LocalSize',[]);	
	parse(p,varargin{:});
	
	%% Check whether the file is intended to be mesh or not
	[~,~,fext] = fileparts(filepath);
	export_geo= strcmpi(fext,'.geo');
	
	%% Locate Gmsh (if meshing is requested)
	if export_geo
		path_to_geo = filepath;						% Dump the geometry at the requested filepath
		wb_title=sprintf('Writing %s',filepath);	% Waitbar title	
	else
		path_to_geo = tempname;						% Dump the geometry to a temporary file
		wb_title='Writing the temporary file...';
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
	end
	
	%% Preprocess geometry if periodicity is requested
	periodic=p.Results.periodic;
	if any(strcmpi(periodic,{'X','Y','both'}))
		[obj, A_seg]=addSymmetry(obj,periodic);
	elseif strcmpi(periodic,'None')
		A_seg=[];
	else
		error('Wrong argument for periodic option. It can be ''x'', ''y'', ''both'' or ''None''.');
	end
	
	
	%% Format data
	elem_type= p.Results.ElementType;
	valid_elem_type={'Wedge','Hex','Tet','Tri','Quad','HexOnly','QuadOnly'};
	if strcmpi(elem_type,'Tetrahedron')
		elem_type='Tet';
	elseif strcmpi(elem_type,'Hexahedron') || strcmpi(elem_type,'Brick')	% 'Brick' is valid for backward compatibility
		elem_type='Hex';
	elseif strcmpi(elem_type,'Triangular')
		elem_type='Tri';
	elseif strcmpi(elem_type,'Quadrangular')
		elem_type='Quad';
	elseif strcmpi(elem_type,'HexahedronOnly') || strcmpi(elem_type,'BrickOnly')
		elem_type='HexOnly';
	elseif strcmpi(elem_type,'QuadrangularOnly')
		elem_type='QuadOnly';
	elseif ~any(strcmpi(elem_type,valid_elem_type)) 
		error('Unrecognized element type. It can be: ''%s''.',strjoin(valid_elem_type,''', '''))
	end
	
	mesh3D = strcmpi(elem_type, 'Hex') || strcmpi(elem_type, 'Wedge')  || strcmpi(elem_type, 'Tet') || strcmpi(elem_type, 'HexOnly');	
	defaultElementSize=p.Results.ElementSize;
	medium=~isequal([0 0 0],p.Results.Medium);
	mediumElementSize=p.Results.MediumElementSize;
	if ~medium && mediumElementSize %	Missing option 'medium'
		error('Specify the size of the embedding medium first with option ''medium''.');
	end
	if defaultElementSize==0
		defaultElementSize=obj.Resolution;	%	Use estimate for spatial resolution
	end		
	slope=p.Results.gradient;
	thickness=p.Results.thickness;
	if thickness==0
		thickness=defaultElementSize;
	elseif ~mesh3D
		error('%s element are not compatible with 2D mesh.', elem_type);
	end


	Curv=p.Results.Curvature;
	if Curv~=0
		if slope~=0
			warning('Non constant element sizes at boundaries are inconsistent with the ''gradient'' option. You might get unexeptected results.');
		end
	end
	if numel(thickness)>1
		warning('The thickness must be a scalar value. I''m using the first value.');
		thickness=thickness(1);
	end

	grainPrefix=p.Results.grainPrefix;

	
	segments=obj.Segments;
	vtx=obj.V;				

	%%	The microstructure is embedded in a medium
	if medium
		ROI=obj.size.ROI;
		dx=p.Results.Medium(1);
		dy=p.Results.Medium(2);
		if mesh3D
			dz=p.Results.Medium(3);
			dmin=min([([dx dy]-ROI)/2 dz-thickness]);	%	Track the minimum distance between ROI and boundaries of the medium
			if dmin<0
				error('The dimensions of the medium must be larger than that of the ROI ([%g %g %g]).',ROI,thickness);
			end			
		else
			dz=0;
			dmin=min(([dx dy]-ROI)/2);	%	Track the minimum distance between ROI and boundaries of the medium
			if dmin<0
				error('The dimensions of the medium must be larger than that of the ROI ([%g %g]).',ROI);
			end				
		end
		xmin=min(vtx(:,1));
		xmax=max(vtx(:,1));
		ymin=min(vtx(:,2));
		ymax=max(vtx(:,2));
		Center=[xmax+xmin ymax+ymin]/2;
		C1=Center+[-dx -dy]/2;
		C2=Center+[-dx +dy]/2;
		C3=Center+[+dx +dy]/2;
		C4=Center+[+dx -dy]/2;
		n_vtx=size(vtx,1);
		vtx=[vtx; C1; C2; C3; C4;];
		n_segments=length(segments);
		new_segments=[n_vtx+1 n_vtx+2; n_vtx+2 n_vtx+3; n_vtx+3 n_vtx+4; n_vtx+4 n_vtx+1];
		new_segments=cast(new_segments,'like',segments{1});	%	The news segments must be of the same type as the other ones.
		segments{n_segments+1}=new_segments(1,:)';			%	Add them to the segment list
		segments{n_segments+2}=new_segments(2,:)';
		segments{n_segments+3}=new_segments(3,:)';
		segments{n_segments+4}=new_segments(4,:)';
		if mediumElementSize==0
			q=1.5;	%	Geometric scale for element size in the medium
			n=log(1+dmin/defaultElementSize*(q-1))/log(q);	%	number of elements with geometric increasing size
			mediumElementSize=defaultElementSize*q^n;
		end						
	end

	%%	Numbering the Line Loops
	[LineLoops,PlaneSurface]=uniqueLoops(obj.Grains);	%	Remove duplicates in Line loops

	%%	Waitbar
	n_vtx=size(vtx,1);
	vtxUsed=ismember(1:n_vtx,vertcat(segments{:}));

	n_segments=length(segments);			
	n_loops=length(LineLoops);			
	n_surfaces=height(obj.Grains);
	n_steps=n_vtx+n_segments+n_surfaces+n_loops;
	step=0;
	set(0,'DefaultTextInterpreter','none');
	h = waitbar(0,path_to_geo,'Name',wb_title,'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
	setappdata(h,'canceling',0)

	%% Dumping the geometry
	ffid = fopen(path_to_geo, 'w');
		%	Heading
		fprintf(ffid,'// File generated with MTEX2Gmsh on %s\n\n',datestr(now));

		%	Mesh parameters
		thicknessName='th';
		defaultElementSizeName='e_def';
		mediumThicknessName='th_med';
		mediumElementSizeName='e_med';
		fprintf(ffid,'// Mesh parameters\n');
		fprintf(ffid,'%s=%g;\n',thicknessName,thickness);
		fprintf(ffid,'%s=%g;\n',defaultElementSizeName,defaultElementSize);
		if medium
			fprintf(ffid,'%s=%g;\n',mediumThicknessName,dz);
			fprintf(ffid,'%s=%g;\n',mediumElementSizeName,mediumElementSize);
			n_steps=n_steps+1;					
		end
		local_size=p.Results.LocalSize;
		n_local_size=size(local_size,1);
		A=false(n_vtx,n_local_size);		
		if n_local_size
			local_size_name=cell(0,1);
			for i=1:n_local_size
				id_grain=local_size(i,1);
				local_size_name_i=sprintf('e_Grain_%i',id_grain);
				local_size_name{i}=local_size_name_i;
				fprintf(ffid,'%s=%g;\n',local_size_name_i,local_size(i,2));
				Out=obj.Grains{id_grain,3};
				Out_segIDs=abs(vertcat(Out{:})); %	Concatenate loop-wise
				In=obj.Grains{id_grain,4};
				In=vertcat(In{:});				%	Concatenate grain-wise
				In_segIDs=abs(vertcat(In{:}));	%	Concatenate loop-wise
				Vtx_ids=obj.Segments([Out_segIDs; In_segIDs]);
				A(unique(vertcat(Vtx_ids{:})),i)=true;
			end
		end

		%	Set Kernel Geometry
		if Curv~=0
			fprintf(ffid,'\nSetFactory("OpenCASCADE");\t // Faster computation of the local curvature\n');
		else
			fprintf(ffid,'\nSetFactory("Built-in");\t // Supports squared BSplines\n');					
		end				

		%	Vertices
		fprintf(ffid,'\n// Vertices\n');		
		for i=1:n_vtx
			step=step+1;
			waitbar(step/n_steps,h,'Vertices coordinates')
			if getappdata(h,'canceling')
				delete(h)
				return
			end
			if vtxUsed(i)
				if medium
					fprintf(ffid,'Point(%i)={%g,%g,0,%s};\n',i,vtx(i,1),vtx(i,2),mediumElementSizeName);	%	If the medium is requested, use the related element size by default. Will be overwritten hereafter.
				else
					if isempty(A) || ~any(A(i,:))
						fprintf(ffid,'Point(%i)={%g,%g,0,%s};\n',i,vtx(i,1),vtx(i,2),defaultElementSizeName);
					else
						grain_ids=find(A(i,:));
						fprintf(ffid,'Point(%i)={%g,%g,0,%s};\n',i,vtx(i,1),vtx(i,2),printMin(local_size_name(grain_ids)));
					end
				end
			end
		end

		%	(B-)Splines
		fprintf(ffid,'\n// Grain boundaries\n');		
		for i=1:n_segments
			step=step+1;
			waitbar(step/n_steps,h,'Sections of boundaries')
			if getappdata(h,'canceling')
				delete(h)
				return
			end
			if length(segments{i})==2
				writeSequence(ffid,'Line',i,segments{i});
			elseif length(segments{i})==3 && Curv~=0
				writeSequence(ffid,'Spline',i,segments{i});	%	OpenCASCADE does not support squared BSpline
			else
				writeSequence(ffid,'BSpline',i,segments{i});				
			end
		end

		%	Loops
		fprintf(ffid,'\n// Closed Loops\n');
		for i=1:n_loops
			step=step+1;
			waitbar(step/n_steps,h,'Line loops')
			if getappdata(h,'canceling')
				delete(h)
				return
			end					
			writeSequence(ffid,'Line Loop',i,LineLoops{i});
		end

		%	Surfaces
		fprintf(ffid,'\n// Grains\n');
		for i=1:n_surfaces
			step=step+1;
			waitbar(step/n_steps,h,'Individual grains')
			if getappdata(h,'canceling')
				delete(h)
				return
			end				    
			writeSequence(ffid,'Plane Surface',i,PlaneSurface{i});
		end

		%	Use quandrangular elements for 2D meshing
		if strcmpi(elem_type, 'Hex') || strcmpi(elem_type, 'Quad') || strcmpi(elem_type, 'HexOnly') || strcmpi(elem_type, 'QuadOnly')
			fprintf(ffid,'\n// Quadrangular elements\n');
			fprintf(ffid,'Recombine Surface{1:%i};\n',n_surfaces);
		end

		%	Extrusions
		if mesh3D
			fprintf(ffid,'\n// 3D geometry\n');
			fprintf(ffid,'Extrude {0,0,%s}{\n\t',thicknessName);
			fprintf(ffid,'Surface{1:%i};\n',n_surfaces);
			fprintf(ffid,'\tLayers{1};');
			if ~strcmpi(elem_type, 'Tet')
				fprintf(ffid,'Recombine;');
			end
			fprintf(ffid,'\n}\n');
		end

		%	Add surrounding medium (if requested)
		if medium
			n_surfaces_tot=n_surfaces+1; %	At least, one more surface exists (surrounding the ROI)
			waitbar(step/n_steps,h,'Adding surrounding medium...')
			fprintf(ffid,'\n// Add surrounding medium\n');

			%	Create a volume beneath the grains
			BL=borderLoop(obj);
			if dz>thickness		%	Add medium below the ROI
				fprintf(ffid,'L[]=Extrude{0,0,%s-%s}{\n\t',thicknessName,mediumThicknessName);	%	Extrude each segment of the outer boundary and keep indices of the resulting segments
				writeSequence(ffid,'Line',[],abs(BL));
				fprintf(ffid,'};\n');
				n_loops=n_loops+1;
				fprintf(ffid,'Line Loop(%i)={',n_loops);			%	Line loop on the opposite side of the ROI
				for i=1:length(BL)
					j=(i-1)*4;	%	Index of the opposite segment from segment i
					if BL(i)>0
						fprintf(ffid,'L[%i]',j);
					else
						fprintf(ffid,'-L[%i]',j);	%	Segments are oriented with respect to their parents
					end
					if i~=length(BL)
						fprintf(ffid,',');
					else
						fprintf(ffid,'};\n');
					end
					if mod(i,50)==0
						fprintf(ffid,'\n\t');
					end
				end
				id_SurfaceLoop=n_surfaces+1;
				writeSequence(ffid,'Plane Surface',id_SurfaceLoop,n_loops);
				fprintf(ffid,'Surface Loop(%i)={\n\t',id_SurfaceLoop);
				for i=1:length(BL)	%	List of side surfaces (results from extrusions)
					j=i*4-3;
					fprintf(ffid,'L[%i],',j);
					if mod(i,50)==0
						fprintf(ffid,'\n\t');
					end
				end
				fprintf(ffid,'\n\t');
				fprintf(ffid,'1:%i\n};\n',n_surfaces+1); %	List of upper faces (original grains)
				writeSequence(ffid,'Volume',n_surfaces+1,id_SurfaceLoop);
				n_surfaces_tot=n_surfaces_tot+2;	%	Two more surfaces are necessaray 
			end

			%	Create a volume surrounding the ROI
			n_loops=n_loops+1;
			writeSequence(ffid,'Line Loop',n_loops,n_segments-3:n_segments);		%	Outer boundaries of the medium					
			n_loops=n_loops+1;
			writeSequence(ffid,'Line Loop',n_loops,BL);							%	Outer boundaries of the ROI/Inner boundaries of the medium
			if mesh3D
				writeSequence(ffid,'Plane Surface',n_surfaces+2,[n_loops-1 n_loops]);%	Upper surface of the medium
			else
				writeSequence(ffid,'Plane Surface',n_surfaces+1,[n_loops-1 n_loops]);%	Only one extra surface in 2D
			end
			if dz>thickness	
				fprintf(ffid,'Extrude {0,0,%s-%s}{\n',thicknessName,mediumThicknessName);
				fprintf(ffid,'\tSurface{%i};\n}\n',n_surfaces+2);
			end
			if mesh3D
				fprintf(ffid,'Extrude {0,0,%s}{\n',thicknessName);
				fprintf(ffid,'\tSurface{%i};\n',n_surfaces+2);
				fprintf(ffid,'\tLayers{1}; Recombine;\n}\n');
			end

			%	Set the correct element size in the ROI
			fprintf(ffid,'Characteristic Length {1:%i} = %s;\n',n_vtx-4,defaultElementSizeName);
		end

		%	Physical groups
		if mesh3D
			groupname='Volume';
		else
			groupname='Surface';
		end		
		Ids=obj.Grains.GrainID(:);
		fprintf(ffid,'\n// Sets\n');
		if all(Ids==(1:n_surfaces)')	%	Grains are numbered subsequently
			waitbar(step/n_steps,h,'Physical volumes');
			fprintf(ffid,'For k In {1:%i}\n',n_surfaces);
			if isempty(grainPrefix)
				fprintf(ffid,'\tPhysical %s(k)={k};\n',groupname);
			elseif isa(grainPrefix,'char')
				fprintf(ffid,'\tPhysical %s(Sprintf("%s%%g",k))={k};\n',groupname,grainPrefix);
			else
				error('grainPrefix must be a string or empty')
			end
			fprintf(ffid,'EndFor\n');
		else							%	Instead, use the ID given by MTEX
			for i=1:n_surfaces
				step=step+1;
				waitbar(step/n_steps,h,'Physical volumes')
				if getappdata(h,'canceling')
					return
				end
				if isempty(grainPrefix)
					fprintf(ffid,'Physical Volume(%i)={%i};\n',grainPrefix,Ids(i),i);
				elseif isa(grainPrefix,'char')
					fprintf(ffid,'Physical Volume("%s_%i")={%i};\n',grainPrefix,Ids(i),i);
				else
					error('grainPrefix must be a string or empty')
				end
			end
		end
		if medium
			if mesh3D
				writeSequence(ffid,sprintf('Physical %s',groupname),'"Medium"',n_surfaces+1:n_surfaces_tot);
			else
				writeSequence(ffid,sprintf('Physical %s',groupname),'"Medium"',n_surfaces+1);
			end
		end
		
		if ~isempty(A_seg)
			fprintf(ffid,'\n// Periodicity conditions\n');
			for i=1:size(A_seg,1)
				fprintf(ffid,'Periodic Line {%i} = {%i};\n',A_seg(i,1),A_seg(i,2));
			end
		end

		%	Mesh 2D
		fprintf(ffid,'\n// Mesh\n');
		fprintf(ffid,'Mesh.SubdivisionAlgorithm=0;\n');		% Turn off subdivision
		if strcmpi(elem_type, 'HexOnly')  || strcmpi(elem_type, 'QuadOnly')
			fprintf(ffid,'Mesh.Algorithm=8;\n');			% Use 'Frontal-Delaunay for quads'
		else
			fprintf(ffid,'Mesh.Algorithm=6;\n');			% Use 'Frontal-Delaunay'
		end		
		fprintf(ffid,'Mesh.CharacteristicLengthExtendFromBoundary=1;\n');
		fprintf(ffid,'Mesh.ElementOrder=%i;\n',p.Results.ElementOrder);
		if Curv~=0
			fprintf(ffid,'Mesh.CharacteristicLengthFromCurvature = 1;\n');
			fprintf(ffid,'Mesh.MinimumCirclePoints = %i; // points per 2*pi\n',Curv);
		end
		if slope~=0				
			fprintf(ffid,'Field[1] = Attractor;\n');
			fprintf(ffid,'Field[1].EdgesList ={1:%i};\n',n_segments);
			fprintf(ffid,'Field[2] = MathEval;\n');
			fprintf(ffid,'Field[2].F = "F1*%g+%g";\n',slope,defaultElementSize(1));
			fprintf(ffid,'Background Field=2;\n');
			fprintf(ffid,'Mesh.CharacteristicLengthExtendFromBoundary=0;\n');
			fprintf(ffid,'Mesh 2;\n');
			fprintf(ffid,'Mesh.CharacteristicLengthExtendFromBoundary=1;\n');				
		end
		if strcmpi(elem_type, 'HexOnly')  || strcmpi(elem_type, 'QuadOnly')
			fprintf(ffid,'Mesh.RecombineAll = 0;\n');
			fprintf(ffid,'Mesh.SaveParametric = 0;\n');
			fprintf(ffid,'Mesh.RecombinationAlgorithm = 0;\n');	% Simple algorithm
			fprintf(ffid,'Mesh.SecondOrderLinear = 1;\n');		% 'we don't have the parametrization of the surface' dixit Geuzaine
			fprintf(ffid,'Mesh.SubdivisionAlgorithm=1;\n');		% Full quad
			fprintf(ffid,'Mesh.MeshSizeFactor=%g;\n', 1);
			fprintf(ffid,'RefineMesh;\n');
			if slope==0.0
				fprintf(ffid,'Mesh 2;\n');
			end
		end		
		if medium
			fprintf(ffid,'Mesh.Algorithm3D=4;\n'); %	Use 'Frontal' algorithm for hybrid structured/unstructured grids
		end
		
		% Mesh 3D
		if mesh3D
			fprintf(ffid,'Mesh 3;\n');
		else
			fprintf(ffid,'Mesh 2;\n');
		end
			
		% Partition mesh
		if p.Results.partition
			fprintf(ffid,'\n// Partition mesh\n');
			fprintf(ffid,'PartitionMesh %i;\n',p.Results.partition);
		end

	fclose(ffid);
	delete(h);

	
	%% Open the file with Gmsh and save mesh
	v=p.Results.verbosity;
	if ~export_geo
		str=sprintf('"%s" "%s" -o "%s" -v %i -save',path_to_gmsh,path_to_geo,filepath,v);
		system(str);
		delete(path_to_geo)	%	delete temp file
	end
	
	if nargout==1
		fh=filepath;
	end
end