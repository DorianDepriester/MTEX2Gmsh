classdef gmshGeo
%	gmshGeo Class for storing descriptions of grains. This class is
%	intended to convert the description given by MTEX into a Gmsh-readable
%	format.
%
%	Single indexing helps to navigate within the aforementioned 
%	description. For instance: obj(5) will select the data of the 5th grain 
%	only. One can also select a series of grains from a given phase. E.g.: 
%	obj('Forsterite') will keep only the data related to the phase named 
%	Forsterite.

    properties
		V=[];				%	Vertices (BSpline knots)
		Segments=cell(0,1); %	Lists of knots, defining the BSplines
		Grains=table;       %	Table summurizing the properties of each grain
		SingularPoints=[];  %	List of singular points (Triple junctions, corners etc.)
		Interfaces=table;	%	Phase-to-phase interfaces
		Resolution=0.0;		%	Approximative spatial resolution
	end
    
    methods
		function G=gmshGeo(grains, varargin)
		%GMSHGEO Object constructor.
		%
		%	GMSHGEO(GRAINS) constructs an instance of class GMSHGEO, from
		%	the object GRAINS (of class grain2d).
		%
		%	The returned object contains the full descriptions of both the
		%	geometries and crystallographic properties of each grain (phase
		%	and orientations).
		%
		%	By default, the Region of Interest (RoI) is supposed to be a
		%	rectangle. Thus, the resulting geometry perfectly fits in a
		%	rectangle. Use
		%		GMSHGEO(GRAINS,'rectangularROI',false)
		%	to drop this feature. 
		%
		%	See also calcGrains, mesh, exportGrainProps.
			if ~isa(grains,'grain2d')
				error('Input argument must be of class grain2d');
			end
			p = inputParser;
			addOptional(p,'rectangularROI',true);
			parse(p,varargin{:});

			[Segmts,OuterLoop,InnerLoops,G.SingularPoints]=computeSegments(grains);
			G.V=grains.V;
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
			G.Resolution=mean(grains.perimeter./cellfun('length',grains.poly));	% Mean vertex-to-vertex distance;
			
			if p.Results.rectangularROI
				G=fixRectangularROI(G);
			end
		end
	end
	
	
	methods (Hidden=true)
		sref=subsref(obj,s)
		
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