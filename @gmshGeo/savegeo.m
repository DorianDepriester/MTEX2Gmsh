function fh=savegeo(obj,filepath,varargin)
%SAVEGEO Save the geometry as an input file for Gmsh (*.geo).
%	SAVEGEO(Object,filepath) saves the geometry in the
%	corresponding file path. The element size for meshing is equal
%	to the EBSD resolution.
%
%	SAVEGEO(...,'ElementSize',s) results in element sizes equal to 
%	s.
%
%	SAVEGEO(...,'Thickness',t) sets an extrusion thickness
%	equal to t (equal to element size by default).
%
%	SAVEGEO(...,'gradient',k) results in elements with size equal
%	to s+k*d (d being the distance from the nearest boundary and s 
%	the default element size).
%
%	SAVEGEO(...,'ElementType',type) sets the element type used
%	for meshing. It can be:
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
%	SAVEGEO(...,'ElementOrder',order) sets the element order. The
%	default value is 1 (i.e. linear elements).
%
%	SAVEGEO(...,'Curvature',np) sets the element sizes to be
%	computed depending on the local curvature (np nodes per 2 pi).
%	np==0 disables this option (default).
%	
%	SAVEGEO(...,'grainPrefix',str) defines the name for the element
%	sets corresponding to grains (Physical Volumes in Gmsh). E.g 
%	SAVEGEO(...,'grainPrefix','grain_') will create volumes named
%	'grain_1', 'grain_2' etc.
%	If the argument is empty, no prefix is given and the physical
%	volumes are just numbered as the grains.
%
%	SAVEGEO(...,'medium',S) embeds the ROI inside a cuboid of size 
%	S=[dx dy dz]. The element size in the medium is	increasing with
%	increasing distance from the ROI. The mesh in the medium is 
%	composed of tetrahedron elements.
%
%	SAVEGEO(...,'medium',S,'mediumElementSize',value) sets the
%	element size at the corners of the medium to the given value.
%
%	h=SAVEGEO(...) returns the full filepath where the geometry has
%	been saved.
%
%	See also mesh.

	version='2.3';	%	MTEX2Gmsh version

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
	parse(p,varargin{:}); 

	defaultElementSize=p.Results.ElementSize;
	medium=~isequal([0 0 0],p.Results.Medium);
	mediumElementSize=p.Results.MediumElementSize;
	if ~medium && mediumElementSize %	Missing option 'medium'
		error('Specify the size of the embedding medium first with option ''medium''.');
	end
	if defaultElementSize==0
		defaultElementSize=obj.evalElementSize;	%	Compute the mean node-to-node distance
	end		
	slope=p.Results.gradient;
	thickness=p.Results.thickness;
	if thickness==0
		thickness=defaultElementSize;
	end

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

	%%	Format file path
	[~,~,fext] = fileparts(filepath);
	if(~strcmpi(fext,'.geo'))
		filepath = [filepath '.geo'];	%	Append the extension if missing
	end	
	segments=obj.Segments;
	vtx=obj.V;				

	%%	The microstructure is embedded in a medium
	if medium
		dx=p.Results.Medium(1);
		dy=p.Results.Medium(2);
		dz=p.Results.Medium(3);
		xmin=min(vtx(:,1));
		xmax=max(vtx(:,1));
		ymin=min(vtx(:,2));
		ymax=max(vtx(:,2));
		ROI=obj.size.ROI;
		dmin=min([([dx dy]-ROI)/2 dz-thickness]);	%	Track the minimum distance between ROI and boundaries of the medium
		if dmin<0
			error('The dimensions of the medium must be larger than that of the ROI ([%g %g %g]).',ROI,thickness);
		end
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
	h = waitbar(0,filepath,'Name','Writing the GEO file...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
	setappdata(h,'canceling',0)

	ffid = fopen(filepath, 'w');
		%%	Heading
		fprintf(ffid,'// File generated with MTEX2Gmsh (v %s) on %s\n\n',version,datestr(now));

		%%	Mesh parameters
		thicknessName='th';
		defaultElementSizeName='e_min';
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

		%%	Set Kernel Geometry
		if Curv~=0
			fprintf(ffid,'\nSetFactory("OpenCASCADE");\t // Faster computation of the local curvature\n');
		else
			fprintf(ffid,'\nSetFactory("Built-in");\t // Supports squared BSplines\n');					
		end				

		%%	Vertices
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
					fprintf(ffid,'Point(%i)={%g,%g,0,%s};\n',i,vtx(i,1),vtx(i,2),defaultElementSizeName);
				end
			end
		end

		%%	(B-)Splines
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

		%%	Loops
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

		%%	Surfaces
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

		%%	Use quandrangular elements for 2D meshing
		if strcmpi(elem_type, 'Hex') || strcmpi(elem_type, 'Quad') || strcmpi(elem_type, 'HexOnly') || strcmpi(elem_type, 'QuadOnly')
			fprintf(ffid,'\n// Quadrangular elements\n');
			fprintf(ffid,'Recombine Surface{1:%i};\n',n_surfaces);
		end

		%%	Extrusions
		if strcmpi(elem_type, 'Hex') || strcmpi(elem_type, 'Wedge')  || strcmpi(elem_type, 'Tet') || strcmpi(elem_type, 'HexOnly')
			fprintf(ffid,'\n// 3D geometry\n');
			fprintf(ffid,'Extrude {0,0,%s}{\n\t',thicknessName);
			fprintf(ffid,'Surface{1:%i};\n',n_surfaces);
			fprintf(ffid,'\tLayers{1};');
			if ~strcmpi(elem_type, 'Tet')
				fprintf(ffid,'Recombine;');
			end
			fprintf(ffid,'\n}\n');
		end

		%%	Add surrounding medium (if requested)
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
			writeSequence(ffid,'Plane Surface',n_surfaces+2,[n_loops-1 n_loops]);%	Upper surface of the medium
			if dz>thickness	
				fprintf(ffid,'Extrude {0,0,%s-%s}{\n',thicknessName,mediumThicknessName);
				fprintf(ffid,'\tSurface{%i};\n}\n',n_surfaces+2);
			end
			fprintf(ffid,'Extrude {0,0,%s}{\n',thicknessName);
			fprintf(ffid,'\tSurface{%i};\n',n_surfaces+2);
			fprintf(ffid,'\tLayers{1}; Recombine;\n}\n');

			%	Set the correct element size in the ROI
			fprintf(ffid,'Characteristic Length {1:%i} = %s;\n',n_vtx-4,defaultElementSizeName);
		end

		%%	Physical volumes
		Ids=obj.Grains.GrainID(:);
		fprintf(ffid,'\n// Sets\n');
		if all(Ids==(1:n_surfaces)')	%	Grains are numbered subsequently
			waitbar(step/n_steps,h,'Physical volumes');
			fprintf(ffid,'For k In {1:%i}\n',n_surfaces);
			if isempty(grainPrefix)
				fprintf(ffid,'\tPhysical Volume(k)={k};\n');
			elseif isa(grainPrefix,'char')
				fprintf(ffid,'\tPhysical Volume(Sprintf("%s%%g",k))={k};\n',grainPrefix);
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
			writeSequence(ffid,'Physical Volume','"Medium"',n_surfaces+1:n_surfaces_tot);
		end

		%%	Mesh
		fprintf(ffid,'\n// Mesh\n');
		fprintf(ffid,'Mesh.CharacteristicLengthExtendFromBoundary=1;\n');
		fprintf(ffid,'Mesh.ElementOrder=%i;\n',p.Results.ElementOrder);
		if strcmpi(elem_type, 'HexOnly')  || strcmpi(elem_type, 'QuadOnly')
			fprintf(ffid,'Mesh.MeshSizeFactor=%d;\n', 0.5);	% A 2.5 factor before subviding keeps the number of nodes almost constant
			fprintf(ffid,'Mesh.SubdivisionAlgorithm=1;\n');
		end
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
		if medium
			fprintf(ffid,'Mesh.Algorithm3D=4;\n'); %	Use 'Frontal' algorithm for hybrid structured/unstructured grids
		end

	fclose(ffid);
	delete(h);

	if nargout==1
		fh=filepath;
	end
end
		
function writeSequence(ffid,title,idx,Seq)
    if ~isempty(Seq)
		if isempty(idx)
		    fprintf(ffid,'%s{%i',title,Seq(1));
		else
			if isnumeric(idx)
				fprintf(ffid,'%s(%i)={%i',title,idx,Seq(1));
			else
				fprintf(ffid,'%s(%s)={%i',title,idx,Seq(1));
			end
		end
		n=length(Seq);
		if n==1
			fprintf(ffid,'};\n');
		else
			pending=[];
			n_pend=0;
			for k=2:n
				if k~=n
					if (Seq(k)==Seq(k-1)+1 || Seq(k)==Seq(k-1)-1) 
						pending=Seq(k);
						n_pend=n_pend+1;
					else
						if isempty(pending)
							fprintf(ffid,',%i',Seq(k));
						elseif n_pend>1
							fprintf(ffid,':%i,%i',pending,Seq(k));
						else
							fprintf(ffid,',%i,%i',pending,Seq(k));
						end
						pending=[];
						n_pend=0;
					end
				else
					if (Seq(k)==Seq(k-1)+1 || Seq(k)==Seq(k-1)-1)
						fprintf(ffid,':%i};\n',Seq(k));
					else
						if isempty(pending)						
							fprintf(ffid,',%i};\n',Seq(k));
						else
							if n_pend>1
								fprintf(ffid,':%i,%i};\n',pending,Seq(k));
							else
								fprintf(ffid,',%i,%i};\n',pending,Seq(k));
							end
						end
					end
				end
			end
		end
    end
end

function segList = borderLoop(G)
	Border=G.Interfaces{'ROI Border',1};
	Border=Border{1};
	F=vertcat(G.Segments{Border});
	F=[F(1:2:end) F(2:2:end)];
	dataType='single';
	p=EulerPath(F,dataType);
	p=p{1};
	nseg=length(p)-1;
	segList=zeros(nseg,1,dataType);			%	Cast the segment list like that of other loops
	for i=1:nseg
		I=find(F(:,1)==p(i) & F(:,2)==p(i+1),1,'first');
		if isempty(I)
			I=find(F(:,2)==p(i) & F(:,1)==p(i+1),1,'first');
			segList(i)=-cast(Border(I),dataType);	%	Border(I) is unsigned 
		else
			segList(i)=Border(I);
		end
	end
	segList=cast(segList,'like',G.Grains.OuterLoop{1});
end

function [LineLoops,PlaneSurface]=uniqueLoops(Grains)
	h=waitbar(0,'Removing loop duplicates','Name','Closed loops');
	n_grains=height(Grains);
	PlaneSurface=cell(height(Grains),1);
	nLoops=n_grains+sum(cellfun(@length,Grains.InnerLoops));	%	Overall number of loops
	LineLoops=cell(nLoops,1);
	jmax=0;	%	Number of unique loops
	for i=1:n_grains
		waitbar(i/n_grains,h);
		OuterLoop=Grains.OuterLoop{i};
		InnerLoops=Grains.InnerLoops{i};
		Seq=zeros(1+length(InnerLoops),1,'int32');
		new=true;
		for j=1:size(LineLoops,1)
			if isempty(LineLoops{j})
				break
			end
			if isequal(sort(abs(OuterLoop)),sort(abs(LineLoops{j})))
				new=false;
				break
			end
		end
		if new
			LineLoops{j}=OuterLoop;
		end
		Seq(1)=j;
		jmax=max(j,jmax);

		for k=1:length(InnerLoops)
			new=true;			
			for j=1:size(LineLoops,1)
				if isempty(LineLoops{j})
					break
				end
				if isequal(sort(abs(InnerLoops{k})),sort(abs(LineLoops{j})))
					new=false;
					break
				end
			end
			if new
				LineLoops{j}=InnerLoops{k};
			end
			Seq(1+k)=j;
			jmax=max(j,jmax);			
		end
		PlaneSurface{i}=Seq;
	end
	LineLoops=LineLoops(1:jmax);
	close(h);
end