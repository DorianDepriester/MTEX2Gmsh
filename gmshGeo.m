classdef gmshGeo
    
    properties
        V=[];
        Segments=cell(0,1);
        Grains=table;
        SingularPoints=[];
		Interfaces=struct;
    end
    
    methods
        function G=gmshGeo(grains)
        %Object constructor.
        % Computes the gmshGeo object from a set of grains (grain2d object)
			if ~isa(grains,'grain2d')
				error('Input argument must be of class grain2d');
			end
            [Segmts,OuterLoop,InnerLoops,G.V,G.SingularPoints]=computeSegments(grains);
			Id=grains.id;
			phaseList=grains.mineralList;
			Phase=phaseList(grains.phaseId)';
			ng=length(grains);
			phi1=zeros(ng,1);Phi=phi1;phi2=phi1;
			convention='Bunge';
			h=waitbar(0,sprintf('Computing Euler angles (with %s convention)',convention),'Name','Individual grain properties');
			for i=1:ng
				waitbar(i/ng,h);
				[phi1(i),Phi(i),phi2(i)]=Euler(grains(i).meanOrientation,convention);
			end
			waitbar(1,h,'Tabular formating');
			G.Grains=table(Id,Phase,OuterLoop,InnerLoops,phi1,Phi,phi2);
			close(h);
			
			
			
			% Initialize the list of phase-to-phase interfaces
			np=length(grains.mineralList);			
			for i=1:np
				for j=i:np
					p1=grains.mineralList{i};
					p2=grains.mineralList{j};
					if ~strcmpi(p1,'notIndexed') || ~strcmpi(p2,'notIndexed')
						strname=genvarname([p1 '_' p2]);
						GB.(strname)=uint16([]);
					end
				end
			end
			border='Border';	% Name for the domain boundary
			GB.(border)=uint16([]);
			
			% Fetch the interfaces
			for ids=1:size(Segmts,1)
				ph2ph=Segmts{ids,2};
				i=min(ph2ph);
				j=max(ph2ph);
				if i==0
					strname=border;
				else
					p1=grains.mineralList{i};
					p2=grains.mineralList{j};					
					strname=genvarname([p1 '_' p2]);
				end
				GB.(strname)=[GB.(strname) ids];
			end
			
			% Remove empty sets
			fn = fieldnames(GB);
			tf = cellfun(@(c) isempty(GB.(c)), fn);
			
			% Update properties
			G.Interfaces = rmfield(GB, fn(tf));
			G.Segments=Segmts(:,1);
        end
        
        function plot(obj,varargin)
        %Plots the segments found in each grains.
        % PLOT(Gmsh) plots all the segments found.
 			if nargin>1
 				GrainIDs=varargin{1};
				lineList=uint16(segmentList(obj,GrainIDs));
			else
				lineList=[];
			end
            pendingLines=true(length(obj.Segments),1);
			intnames=fieldnames(obj.Interfaces);
			nInt=length(intnames);
			C=lines(nInt);
			for i=1:nInt
				intname=intnames{i};
				lineSet=obj.Interfaces.(intname);
				X=[];Y=[];
				for j=1:length(lineSet)
					lineID=lineSet(j);
					if isempty(lineList) || ismember(lineID,lineList)
						pendingLines(lineID)=0;
						Vids=obj.Segments{lineID,1};
						x=obj.V(Vids,1);
						y=obj.V(Vids,2);
						if Vids(1)==Vids(end)
							XYbs=BSpline([x(1:end-1),y(1:end-1)],'order',2,'periodic',true);
						else
							XYbs=BSpline([x,y],'order',2);
						end
						X=[X; NaN; XYbs(:,1)];
						Y=[Y; NaN; XYbs(:,2)];
					end
				end
				p=plot(X,Y);
				set(p,'Color',C(i,:),'LineWidth',2);
				hold on
			end
            hold off
            axis equal
			h=legend(intnames,'Interpreter', 'none');
 			fontSize=14;
 			set(h,'FontSize',fontSize);
			xlabel('e_1','Interpreter', 'tex');
			ylabel('e_2','Interpreter', 'tex');
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
        end
        
        function savegeo(obj,filepath,varargin)
        %Saves the geometry as an input file for Gmsh (*.geo).
        %	- SAVEGEO(Gmsh,filepath) saves the geometry in the
        % corresponding file path. The element size for meshing is equalt
        % to the EBSD resolution.
		%	- SAVEGEO(...,'ElementSize',value) results in element sizes
		%  equal to the given value.
		%	- SAVEGEO(...,'Thickness',value) sets an extrusion thickness
		%  equal to the given value (equal to element size by default).
        %	- SAVEGEO(...,'gradient',k) results in elements with size equal
        %  to s+k*d (d being the distance from the nearest boundary and s 
		%  the default element size).
		%	- SAVEGEO(...,'ElementType',type) sets the element type used
		%  for meshing. It can be:
		%		-'Wedge' (default) for Wedge elements,
		%		-'Brick' for quadrangular (2D)/Brick (3D) elements.
		%	- SAVEGEO(...,'Curvature',np) sets the element sizes to be
		%	computed depending on the local curvature (np nodes per 2 pi).
		%	np==0 disables this option (default).
		%	- SAVEGEO(...,'medium',S) embeds the ROI inside a
		%	cuboid of size S=[dx dy dz]. The element size in the medium is
		%	increasing with increasing distance from the ROI. The mesh in
		%	the medium is composed of tetrahedron elements.
		%	- SAVEGEO(...,'medium',S,'mediumElementSize',value) sets the
		%   element size at the corners of the medium to the given value.
			version='1.0';	% MTEX2Gmsh version
		
			%% Parse optional parameters
            p = inputParser;
            addOptional(p,'ElementSize',0);
            addOptional(p,'thickness',0);
            addOptional(p,'gradient',0);
            addOptional(p,'ElementType','Wedge');
            addOptional(p,'Curvature',0);
            addOptional(p,'Medium',[0 0 0]);
            addOptional(p,'MediumElementSize',0);
            parse(p,varargin{:}); 
			
            defaultElementSize=p.Results.ElementSize;
			medium=~isequal([0 0 0],p.Results.Medium);
			mediumElementSize=p.Results.MediumElementSize;
			if ~medium && mediumElementSize % Missing option 'medium'
				error('Specify the size of the embedding medium first with option ''medium''.');
			end
			if defaultElementSize==0
				defaultElementSize=obj.evalElementSize;	% Compute the mean node-to-node distance
			end		
            slope=p.Results.gradient;
            thickness=p.Results.thickness;
			if thickness==0
				thickness=defaultElementSize;
			end
			elem_type= p.Results.ElementType;
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
			
			%% Format file path
            [~,~,fext] = fileparts(filepath);
			if(~strcmpi(fext,'.geo'))
                filepath = [filepath '.geo'];	% Append the extension if missing
			end	
            segments=obj.Segments;
            vtx=obj.V;				

			%% The microstructure is embedded in a medium
			if medium
				dx=p.Results.Medium(1);
				dy=p.Results.Medium(2);
				dz=p.Results.Medium(3);
				xmin=min(vtx(:,1));
				xmax=max(vtx(:,1));
				ymin=min(vtx(:,2));
				ymax=max(vtx(:,2));
				ROI=obj.size.ROI;
				dmin=min([([dx dy]-ROI)/2 dz-thickness]);	% Track the minimum distance between ROI and boundaries of the medium
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
				new_segments=cast(new_segments,'like',segments{1});	% The news segments must be of the same type as the other ones.
				segments{n_segments+1}=new_segments(1,:)';			% Add them to the segment list
				segments{n_segments+2}=new_segments(2,:)';
				segments{n_segments+3}=new_segments(3,:)';
				segments{n_segments+4}=new_segments(4,:)';
				if mediumElementSize==0
					q=1.5;	% Geometric scale for element size in the medium
					n=log(1+dmin/defaultElementSize*(q-1))/log(q);	% number of elements with geometric increasing size
					mediumElementSize=defaultElementSize*q^n;
				end						
			end
			
			%% Numbering the Line Loops
			[LineLoops,PlaneSurface]=uniqueLoops(obj.Grains);	% Remove duplicates in Line loops
			
            %% Waitbar
            n_vtx=size(vtx,1);
            vtxUsed=ismember(1:n_vtx,cell2mat(segments));
			
            n_segments=length(segments);			
            n_loops=length(LineLoops);			
            n_surfaces=height(obj.Grains);
            n_steps=n_vtx+n_segments+n_surfaces+n_loops;
            step=0;
            set(0,'DefaultTextInterpreter','none');
            h = waitbar(0,filepath,'Name','Writing the GEO file...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            setappdata(h,'canceling',0)

            ffid = fopen(filepath, 'w');
				%% Heading
				fprintf(ffid,'// File generated with MTEX2Gmsh (v %s) on %s\n\n',version,datestr(now));
				
                %% Mesh parameters
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
				
				%% Set Kernel Geometry
				if Curv~=0
					fprintf(ffid,'\nSetFactory("OpenCASCADE");\t // Faster computation of the local curvature\n');
				else
					fprintf(ffid,'\nSetFactory("Built-in");\t // Supports squared BSplines\n');					
				end				

                %% Vertices
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
							fprintf(ffid,'Point(%i)={%g,%g,0,%s};\n',i,vtx(i,1),vtx(i,2),mediumElementSizeName);	% If the medium is requested, use the related element size by default. Will be overwritten hereafter.
						else
							fprintf(ffid,'Point(%i)={%g,%g,0,%s};\n',i,vtx(i,1),vtx(i,2),defaultElementSizeName);
						end
                    end
                end

                %% (B-)Splines
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
						writeSequence(ffid,'Spline',i,segments{i});	% OpenCASCADE does not support squared BSpline
					else
						writeSequence(ffid,'BSpline',i,segments{i});				
					end
				end
				
				%% Loops
				fprintf(ffid,'\n// Closed Loops\n');
				for i=1:n_loops
					step=step+1;
					waitbar(step/n_steps,h,'Line loops')
					if getappdata(h,'canceling')
                        delete(h)
                        return
					end					
					writeSequence(ffid,'Line Loop',10*i,LineLoops{i});	% label times 10 in order to avoid label conflict with OpenCASCADE (bug)
				end

                %% Surfaces
                fprintf(ffid,'\n// Grains\n');
				for i=1:n_surfaces
                    step=step+1;
                    waitbar(step/n_steps,h,'Individual grains')
                    if getappdata(h,'canceling')
                        delete(h)
                        return
                    end                    
					writeSequence(ffid,'Plane Surface',i,PlaneSurface{i}*10);	% Times 10 to keep consistent with with previous hack
				end
                
                %% Use quandrangular elements for 2D meshing
                if strcmpi(elem_type, 'Brick')
                    fprintf(ffid,'\n// Quadrangular elements\n');
					fprintf(ffid,'Recombine Surface{1:%i};\n',n_surfaces);
                end

                %% Extrusions
				fprintf(ffid,'\n// 3D geometry\n');
				fprintf(ffid,'Extrude {0,0,%s}{\n\t',thicknessName);
				fprintf(ffid,'Surface{1:%i};\n',n_surfaces);
				fprintf(ffid,'\tLayers{1};');
				if ~strcmpi(elem_type, 'Tet') && ~strcmpi(elem_type, 'Tetrahedron')
					fprintf(ffid,'Recombine;');
				end
				fprintf(ffid,'\n}\n');
				
				%% Add surrounding medium (if requested)
				if medium
					n_surfaces_tot=n_surfaces+1; % At least, one more surface exists (surrounding the ROI)
					waitbar(step/n_steps,h,'Adding surrounding medium...')
					fprintf(ffid,'\n// Add surrounding medium\n');
					
					% Create a volume beneath the grains
					BL=borderLoop(obj);
					if dz>thickness		% Add medium below the ROI
						fprintf(ffid,'L[]=Extrude{0,0,%s-%s}{\n\t',thicknessName,mediumThicknessName);	% Extrude each segment of the outer boundary and keep indices of the resulting segments
						writeSequence(ffid,'Line',[],abs(BL));
						fprintf(ffid,'};\n');
						n_loops=n_loops+1;
						fprintf(ffid,'Line Loop(%i)={',n_loops*10);			% Line loop on the opposite side of the ROI
						for i=1:length(BL)
							j=(i-1)*4;	% Index of the opposite segment from segment i
							if BL(i)>0
								fprintf(ffid,'L[%i]',j);
							else
								fprintf(ffid,'-L[%i]',j);	% Segments are oriented with respect to their parents
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
						writeSequence(ffid,'Plane Surface',id_SurfaceLoop,n_loops*10);
						fprintf(ffid,'Surface Loop(%i)={\n\t',id_SurfaceLoop);
						for i=1:length(BL)	% List of side surfaces (results from extrusions)
							j=i*4-3;
							fprintf(ffid,'L[%i],',j);
							if mod(i,50)==0
								fprintf(ffid,'\n\t');
							end
						end
						fprintf(ffid,'\n\t');
						fprintf(ffid,'1:%i\n};\n',n_surfaces+1); % List of upper faces (original grains)
						writeSequence(ffid,'Volume',n_surfaces+1,id_SurfaceLoop);
						n_surfaces_tot=n_surfaces_tot+2;	% Two more surfaces are necessaray 
					end
					
					% Create a volume surrounding the ROI
					n_loops=n_loops+1;
					writeSequence(ffid,'Line Loop',n_loops*10,n_segments-3:n_segments);		% Outer boundaries of the medium					
					n_loops=n_loops+1;
					writeSequence(ffid,'Line Loop',n_loops*10,BL);							% Outer boundaries of the ROI/Inner boundaries of the medium
					writeSequence(ffid,'Plane Surface',n_surfaces+2,[n_loops-1 n_loops]*10);% Upper surface of the medium
					if dz>thickness	
						fprintf(ffid,'Extrude {0,0,%s-%s}{\n',thicknessName,mediumThicknessName);
						fprintf(ffid,'\tSurface{%i};\n}\n',n_surfaces+2);
					end
					fprintf(ffid,'Extrude {0,0,%s}{\n',thicknessName);
					fprintf(ffid,'\tSurface{%i};\n',n_surfaces+2);
					fprintf(ffid,'\tLayers{1}; Recombine;\n}\n');
					
					% Set the correct element size in the ROI
					fprintf(ffid,'Characteristic Length {1:%i} = %s;\n',n_vtx-4,defaultElementSizeName);
				end
				
                %% Physical volumes
				grainPrefix='Grain';
				Ids=obj.Grains.Id(:);
                fprintf(ffid,'\n// Sets\n');
				if all(Ids==(1:n_surfaces)')	% Grains are numbered subsequently
					waitbar(step/n_steps,h,'Physical volumes');
					fprintf(ffid,'For k In {1:%i}\n',n_surfaces);
					fprintf(ffid,'\tPhysical Volume(Sprintf("%s_%%g",k))={k};\n',grainPrefix);
					fprintf(ffid,'EndFor\n');
				else							% Instead, use the ID given by MTEX
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

                %% Mesh
                fprintf(ffid,'\n// Mesh\n');
                fprintf(ffid,'Mesh.CharacteristicLengthExtendFromBoundary=1;\n');
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
					fprintf(ffid,'Mesh.Algorithm3D=4;\n'); % Use 'Frontal' algorithm for hybrid structured/unstructured grids
				end
				
            fclose(ffid);
            delete(h);            
		end
        
        function plotElementSize(obj,minSize,slope,varargin)
        %Plots the map of the element size when gradient is enabled.
        %   -PLOTELEMENTSIZE(obj,minSize,slope) computes the field
        % corresponding to minSize+slope*d, with d the distance from the
        % nearest vertex. Then it plots it as a 2D map.
        %   -PLOTELEMENTSIZE(...,'samples',n) uses n samples in each
        %   directions (default is 200).
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
		end      

		function exportGrainProps(obj,filename)
		%Exports grain properties (IDs, phase and orientation) as ASCII
		%data in a CSV file.
		% EXPORTGRAINPROPS(G,'filename') exports grain properties stored
		% in the G object (of class gmshGeo) in the ASCII file named
		% 'filename'.
			set(0,'DefaultTextInterpreter','none');
			h = waitbar(0,'Grain properties','Name','Writing the CSV file...');
			ffid = fopen(filename, 'w');
			fprintf(ffid,'"GrainID"\t"Phase"\t"phi1"\t"Phi"\t"phi2"\n');
			G=obj.Grains;
			ng=height(G);
			for i=1:ng
				waitbar(i/ng,h);
				Phase=G.Phase(i);
				fprintf(ffid,'%i\t"%s"\t%g\t%g\t%g\n',G.Id(i),Phase{1},G.phi1(i),G.Phi(i),G.phi2(i));
			end
			fclose(ffid);
			close(h);
		end
		
		function s=evalElementSize(obj)
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
		
		function segmts=simplify(obj,varargin)
			%%Applies the Douglas-Peucker algorithm to reduce the number of
			%%points in each segment.
			%   SIMPLIFY(G) reduces the number elements with penalty lenght
			% equal to one tenth of the default element size.
			%   SIMPLIFY(G,epsilon) uses epsilon as the penalty length.
			if nargin==1
				epsilon=obj.evalElementSize/10;
			else
				epsilon=varargin{1};
			end
			segmts=obj.Segments;
			for i=1:length(obj.Segments)
				segmt=obj.Segments{i};
				remains=DouglasPeucker(obj.V(segmt,:),epsilon);
				segmts{i}=segmt(remains);
			end
		end

		function s=size(obj)
			s.numberOfGrains=height(obj.Grains);
			vtx=obj.V;
			xmin=min(vtx(:,1));
			xmax=max(vtx(:,1));
			ymin=min(vtx(:,2));
			ymax=max(vtx(:,2));
			s.ROI=[xmax-xmin ymax-ymin];
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
		if X(1) == X(end) && Y(1)==Y(end)	% Closed loop
			mid=round(npt/2);				% Add mid-point
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

function  segIDs  = segmentList(G,GrainIDs)
	t=ismember(G.Grains.Id,GrainIDs);
	%% Outer loops
	Out=G.Grains{t,3};
	Out_segIDs=abs(vertcat(Out{:})); % Concatenate loop-wise
	%% Inner loops
	In=G.Grains{t,4};
	In=vertcat(In{:});				% Concatenate grain-wise
	In_segIDs=abs(vertcat(In{:}));	% Concatenate loop-wise 
	segIDs=unique([Out_segIDs; In_segIDs]);
end

function segList = borderLoop(G)
	Border=G.Interfaces.Border;
	F=vertcat(G.Segments{Border});
	F=[F(1:2:end) F(2:2:end)];
	dataType='single';
	p=EulerPath(F,dataType);
	p=p{1};
	nseg=length(p)-1;
	segList=zeros(nseg,1,dataType);			% Cast the segment list like that of other loops
	for i=1:nseg
		I=find(F(:,1)==p(i) & F(:,2)==p(i+1),1,'first');
		if isempty(I)
			I=find(F(:,2)==p(i) & F(:,1)==p(i+1),1,'first');
			segList(i)=-cast(Border(I),dataType);	% Border(I) is unsigned 
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
	nLoops=n_grains+sum(cellfun(@length,Grains.InnerLoops));	% Overall number of loops
	LineLoops=cell(nLoops,1);
	jmax=0;	% Number of unique loops
	for i=1:n_grains
		waitbar(i/n_grains,h);
		OuterLoop=Grains.OuterLoop{i};
		InnerLoops=Grains.InnerLoops{i};
		Seq=zeros(1+length(InnerLoops),1,'int32');
		new=true;
		for j=1:size(LineLoops,1);
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
			for j=1:size(LineLoops,1);
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
