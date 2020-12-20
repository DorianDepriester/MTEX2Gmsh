classdef gmshGeo
    
    properties
		V=[];		       %	Vertices (BSpline knots)
		Segments=cell(0,1); %	Lists of knots, defining the BSplines
		Grains=table;       %	Table summurizing the properties of each grain
		SingularPoints=[];  %	List of singular points (Triple junctions, corners etc.)
		Interfaces=table;  %	Phase-to-phase interfaces
	end
    
    methods
		function G=gmshGeo(grains)
		%GMSHGEO Object constructor.
		%
		%	GMSHGEO(GRAINS) constructs an instance of class GMSHGEO, from
		%	the object GRAINS (of class grain2d).
		%
		%	The returned object contains the full descriptions of both the
		%	geometries and crystallographic properties of each grain (phase
		%	and orientations).
		%
		%	Note that single indexing helps to navigate within those
		%	descriptions. For instance: obj(5) will select the data of the
		%	5th grain only. One can also select a series of grains from a
		%	given phase. E.g.: obj('Forsterite') will keep only the data
		%	related to the phase named Forsterite.
		%
		%	See also calcGrains, mesh, exportGrainProps.
      if ~isa(grains,'grain2d')
        error('Input argument must be of class grain2d');
      end
    
      [Segmts,OuterLoop,InnerLoops,G.SingularPoints]=computeSegments(grains);
			G.V=grains.boundary.V;
			GrainID=grains.id;
			phaseList=grains.mineralList;
			Phase=phaseList(full(grains.phaseId))';
						
			convention='Bunge';
      [phi1,Phi,phi2] = Euler(grains.meanRotation,convention);
			
			h = waitbar(1,'Tabular formating');
			G.Grains=table(GrainID,Phase,OuterLoop,InnerLoops,phi1,Phi,phi2);
			close(h);
			
			T=table([],'VariableNames',{'SegmentIDs'});
			for ids=1:size(Segmts,1)
				ph2ph=Segmts{ids,2};
				i=min(ph2ph);
				j=max(ph2ph);
				if i==0
					name='ROI Border';
				else
					p1=grains.mineralList{i};
					p2=grains.mineralList{j};					
					name=[p1 ' - ' p2];
				end
				row=find(strcmp(T.Properties.RowNames,name));
				if isempty(row)
					newRow=table({ids},'VariableNames',{'SegmentIDs'},'RowNames',{name});
					T=[T; newRow]; %#ok<AGROW>
				else
					T.SegmentIDs{row}=cat(1,T.SegmentIDs{row},ids);
				end
			end
			G.Interfaces=T;
			G.Segments=Segmts(:,1);
		end
		
		function plot(obj,varargin)
		%PLOT Plot the segments found in each grains.
		%	PLOT(obj) plots all the segments.
		%
		%	PLOT(obj(I)) plots the segments of grains whose indices are
		%	given by the array I.
		%
		%	PLOT(obj(P)) plots the segments of grains of phase P only,
		%	where P is a string.
		%
		%	See also plotElementSize
			if isempty(obj.Grains) || isequal(obj.Grains,struct)
				warning('Empty set. Nothing to plot.')
				return
			end
			T=obj.Interfaces;
			intnames=T.Properties.RowNames;
			nInt=length(intnames);
			D=cell(2*nInt,1);
			for i=1:nInt
				lineSet=T{i,1};
				lineSet=vertcat(lineSet{:});
				X=[];Y=[];
				for j=1:length(lineSet)
					lineID=lineSet(j);
					Vids=obj.Segments{lineID,1};
					x=obj.V(Vids,1);
					y=obj.V(Vids,2);
					if Vids(1)==Vids(end)
						XYbs=BSpline([x(1:end-1),y(1:end-1)],'order',2,'periodic',true);
					else
						XYbs=BSpline([x,y],'order',2);
					end
					X=[X; NaN; XYbs(:,1)]; %#ok<AGROW>
					Y=[Y; NaN; XYbs(:,2)]; %#ok<AGROW>
				end
				D{2*i-1}=X;
				D{2*i}=Y;
			end
		    p=plot(D{:},varargin{:});
			set(p,'LineWidth',2);
			for i=1:nInt
				set(p(i),'DisplayName',intnames{i})
			end
			xlabel('x');
			ylabel('y');
			setPlotOrientation
		end
		
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
		%		-'Wedge' (default) for 6-node 3D elements,
		%		-'Hex' for 8-node 3D elements,
		%		-'Tet' or 'Tetrahedron' for 4-node 3D elements,
		%		-'Tri' or 'Triangular' for 3-node 2D elements,
		%		-'Quad' or 'Quadrangular' for 4-node 2D elements,
		%		-'HexOnly' for hexahedron elements only (no tet),
		%		-'QuadOnly' for quandrangular elements only (no triangle).
		%
		%	SAVEGEO(...,'ElementOrder',order) sets the element order. The
		%	default value is 1 (i.e. linear elements).
		%
		%	SAVEGEO(...,'Curvature',np) sets the element sizes to be
		%	computed depending on the local curvature (np nodes per 2 pi).
		%	np==0 disables this option (default).
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
		
			version='2.2';	%	MTEX2Gmsh version
		
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
				grainPrefix='Grain';
				Ids=obj.Grains.GrainID(:);
				fprintf(ffid,'\n// Sets\n');
				if all(Ids==(1:n_surfaces)')	%	Grains are numbered subsequently
					waitbar(step/n_steps,h,'Physical volumes');
					fprintf(ffid,'For k In {1:%i}\n',n_surfaces);
					fprintf(ffid,'\tPhysical Volume(Sprintf("%s_%%g",k))={k};\n',grainPrefix);
					fprintf(ffid,'EndFor\n');
				else							%	Instead, use the ID given by MTEX
					for i=1:n_surfaces
						step=step+1;
						waitbar(step/n_steps,h,'Physical volumes')
						if getappdata(h,'canceling')
							return
						end				    
						fprintf(ffid,'Physical Volume("%s_%i")={%i};\n',grainPrefix,Ids(i),i);
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
					fprintf(ffid,'Mesh.MeshSizeFactor=%d;\n', 2.5);	% A 2.5 factor before subviding keeps the number of nodes almost constant
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
		
		function mesh(obj,outputFilePath,varargin)
			%MESH Mesh the geometry and export the mesh into the
			%requested file.
			%
			%	MESH(obj,outputFile) meshes the geometry and saves the 
			%	results at the specified location. If the latter is only a 
			%	filename (no path), the mesh is saved in the current 
			%	working directory. The element size is roughly equal to the
			%	EBSD resolution.
			%
			%	MESH(...,'ElementSize',s) results in element sizes
			%	equal to s.
			%
			%	MESH(...,'Thickness',t) sets an extrusion thickness equal
			%	to T (equal to element size by default).
			%
			%	MESH(...,'gradient',k) results in elements with size equal
			%	to s+k*d (d being the distance from the nearest boundary 
			%	and s the default element size).
			%
			%	MESH(...,'ElementType',type) sets the element type used
			%	for meshing. It can be:
			%		-'Wedge' (default) for 6-node 3D elements,
			%		-'Hex' for 8-node 3D elements,
			%		-'Tet' or 'Tetrahedron' for 4-node 3D elements,
			%		-'Tri' or 'Triangular' for 3-node 2D elements,
			%		-'Quad' or 'Quadrangular' for 4-node 2D elements.
			%
			%	MESH(...,'Curvature',np) sets the element sizes to be
			%	computed depending on the local curvature (np nodes per 2 
			%	pi). np==0 disables this option (default).
			%
			%	MESH(...,'medium',S) embeds the ROI inside a cuboid of size 
			%	S=[dx dy dz]. The element size in the medium is	increasing
			%	with increasing distance from the ROI. The mesh in the 
			%	medium is composed of tetrahedron elements.
			%
			%	MESH(...,'medium',S,'mediumElementSize',value) sets the 
			%	element size at the corners of the medium to the given 
			%	value.
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
		
		function plotElementSize(obj,minSize,slope,varargin)
		%PLOTELEMENTSIZE Plots the map of the element size when gradient is 
		%enabled.
		%
		%	PLOTELEMENTSIZE(obj,minSize,slope) computes the field
		%	corresponding to minSize+slope*d, with d the distance from the
		%	nearest vertex. Then it plots it as a 2D map.
		%
		%	PLOTELEMENTSIZE(...,'samples',n) uses n samples in each
		%	directions (default is 200).
		%
		%	This function is intended to check whether the slope value for
		%	writing the .geo file is correct fits the user's needs.
		%
		%	See also plot, savegeo.
		    p = inputParser;
			addOptional(p,'samples',200);
		    parse(p,varargin{:});
		    npt=p.Results.samples;

		    vtx=obj.V;
		    Xmin=min(vtx(:,1));
		    Xmax=max(vtx(:,1));
		    Ymin=min(vtx(:,2));
		    Ymax=max(vtx(:,2));

		    Xlin=linspace(Xmin,Xmax,npt);
		    Ylin=linspace(Ymin,Ymax,npt);
		    [X,Y]=meshgrid(Xlin,Ylin);
		    dist=inf(npt);
		    h=waitbar(0,'Computing the distances from each vertex...');
		    nV=size(vtx,1);
		    for i=1:nV
				waitbar(i/nV,h);
				disti=(X-vtx(i,1)).^2+(Y-vtx(i,2)).^2;
				dist=min(dist,disti);
		    end
		    dist=sqrt(dist);
		    elemSize=minSize+slope*dist;
		    close(h)
		    imagesc(Xlin,Ylin,elemSize);
		    colorbar
			axis equal
			xlabel('x');
			ylabel('y');
			setPlotOrientation			
		end      

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
				rodr=Rodrigues(rot);
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
		
		function s=evalElementSize(obj)
			%EVALELEMENTSIZE Automatically evaluate the element size.
			%	EVALELEMENTSIZE(Object) computes the mean node-to-node
			%	distance in each segment.
			%
			%	Note: the value from EVALELEMENTSIZE is used by default in
			%	the savegeo method.
			%
			%	See also savegeo.
			segmts=obj.Segments;
			nseg=length(segmts);
			d=zeros(nseg,1);
			wt=zeros(nseg,1);
			for i=1:nseg
				if length(segmts{i})>2
					XY=obj.V(segmts{i},:);
					dl=sum((XY(2:end,:)-XY(1:end-1,:)).^2,2);
					d(i)=mean(sqrt(dl));
					wt(i)=length(dl);
				end
			end
			keep= d~=0;
			d=d(keep);
			wt=wt(keep);
			s=d'*wt/sum(wt);
		end
		
		function G=simplify(obj,varargin)
			%SIMPLIFY Applies the Douglas-Peucker algorithm to reduce the number of
			%points in each segment.
			%
			%	SIMPLIFY(Object) reduces the number elements with penalty 
			%	lenght equal to one tenth of the default element size.
			%
			%	SIMPLIFY(Object,epsilon) uses epsilon as the penalty
			%	length.
			%
			%	See also plot.
			G=obj;
			if nargin==1
				epsilon=obj.evalElementSize/10;
			else
				epsilon=varargin{1};
			end
			
			%% Apply the Douglas-Peucker algorithm on each segment
			segments=obj.Segments;
			h=waitbar(0,'Applying the Douglas-Peucker algorithm','Name','Simplification of the boundaries','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
		    setappdata(h,'canceling',0)
			n_seg=length(segments);
			for i=1:n_seg
				if getappdata(h,'canceling')
					delete(h)
					return
				end
				waitbar(i/n_seg,h);
				segmt=obj.Segments{i};
				remains=DouglasPeucker(obj.V(segmt,:),epsilon);
				segments{i}=segmt(remains);
			end
			
			%% Remove unused vertices and update the segments
			Vtx=G.V;
			used=ismember(1:size(Vtx,1),cell2mat(segments));	% Track wether the vertices are used or not in the segments
			new_idx=cast(cumsum(used),'like',segments{1});
			sp=G.SingularPoints;
			for j=1:length(segments)
				segments{j}=new_idx(segments{j})';
			end
			G.Segments=segments;			% Update the segments
			G.V=Vtx(used,:);				% Remove unused vertices
			G.SingularPoints=new_idx(sp);	% Update singular points
			delete(h)
		end

		function s=size(obj)
		    %SIZE Number of grains and size of the ROI.
			s.numberOfGrains=height(obj.Grains);
			vtx=obj.V;
			xmin=min(vtx(:,1));
			xmax=max(vtx(:,1));
			ymin=min(vtx(:,2));
			ymax=max(vtx(:,2));
			s.ROI=[xmax-xmin ymax-ymin];
		end	
		
	end
	
	methods (Hidden=true)
		function sref=subsref(obj,s)
		   %	obj(i) only selects the data related to the i-th grain
			switch s(1).type
				case '.'
					sref=builtin('subsref',obj,s);
				case '()'
					sref=obj;
					grain_tab=sref.Grains;
					k = s.subs;
					if length(k)>1
						error('Only single index can be used here. Consider using an array of indices instead.')
					end						

					%%	Select the grains in the table
					if all(cellfun(@(x) isnumeric(x),k))
						rows=k{:};
						if any(rows<1)
							error('Indices should be stricly positive.')
						end
						if any(rows>height(grain_tab))
							error('The index is larger than the number of grains (%i here).',height(grain_tab))
						end
					else
						rows=strcmp(grain_tab.Phase,k);
					end
					grain_tab=grain_tab(rows,:);
					sref.Grains=grain_tab;

					%%	Keep only the related segments
					if isempty(grain_tab)
						sref.Interfaces=[];
					else
						Out=grain_tab{:,3};
						Out_segIDs=abs(vertcat(Out{:})); %	Concatenate loop-wise
						In=grain_tab{:,4};
						In=vertcat(In{:});				%	Concatenate grain-wise
						In_segIDs=abs(vertcat(In{:}));	%	Concatenate loop-wise 
						segIDs=unique([Out_segIDs; In_segIDs]);

						%%	Update the interfaces
						interfaces=sref.Interfaces;
						for i=1:height(interfaces)
							interface=cast(interfaces.SegmentIDs{i},'like',segIDs);
							interfaces.SegmentIDs{i}=interface(ismember(interface,segIDs));
						end
						t=cellfun('isempty',interfaces.SegmentIDs);	% Keep only non empty sets of segments
						sref.Interfaces=interfaces(~t,:);
					end
				case '{}'
					error('gmshGeo:subsref',...
						'Not a supported subscripted reference')
			end
		end
		
		function k = end(obj,~,n)
		%%Overload the end method
			if n>1
				error('Only single index can be used.')
			else
				k=height(obj.Grains);
			end
		end
		
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

function remains=DouglasPeucker(V,epsilon)
	npt=size(V,1);
	if npt>2
		X=V(:,1);
		Y=V(:,2);		
		Pdg=zeros(npt,1);
		Pdg(1)=1;
		Pdg(end)=npt;
		if X(1) == X(end) && Y(1)==Y(end)	%	Closed loop
			mid=round(npt/2);				%	Add mid-point
			Pdg(mid)=mid;
		end
		d=inf;
		imax=0;
		while d>epsilon
			d=0;
			for i=2:(npt-1)
				if Pdg(i)==0
					im=find(Pdg<i & Pdg~=0,1,'last');
					ip=find(Pdg>i & Pdg~=0,1,'first');
					U=[X(ip)-X(im)
					  Y(ip)-Y(im)];
					U=U/norm(U);
					A=[X(i)-X(im)
					  Y(i)-Y(im)];
					di=abs(det([U A]));
					if di>d
						d=di;
						imax=i;
					end
				end
			end
			if imax~=0
				Pdg(imax)=imax;
			end
		end
		remains=Pdg~=0;
	else
		remains=true(npt,1);
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

function setPlotOrientation
	if strcmpi(getMTEXpref('xAxisDirection'),'east') && strcmpi(getMTEXpref('zAxisDirection'),'intoPlane')
		set(gca,'Ydir','reverse')
	elseif strcmpi(getMTEXpref('xAxisDirection'),'north')
		if strcmpi(getMTEXpref('zAxisDirection'),'outOfPlane')
			set(gca,'Ydir','reverse')
		end
		view([90 -90])
	elseif strcmpi(getMTEXpref('xAxisDirection'),'west')
		set(gca,'Xdir','reverse')
		if strcmpi(getMTEXpref('zAxisDirection'),'outOfPlane')
			set(gca,'Ydir','reverse')
		end
	elseif strcmpi(getMTEXpref('xAxisDirection'),'south')
		view([90 -90])
		set(gca,'Xdir','reverse')
		if strcmpi(getMTEXpref('zAxisDirection'),'intoPlane')
			set(gca,'Ydir','reverse')
		end
	end
	axis equal
end