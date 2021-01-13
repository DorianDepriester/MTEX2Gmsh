classdef gmshGeo
    
    properties
		V=[];				%	Vertices (BSpline knots)
		Segments=cell(0,1); %	Lists of knots, defining the BSplines
		Grains=table;       %	Table summurizing the properties of each grain
		SingularPoints=[];  %	List of singular points (Triple junctions, corners etc.)
		Interfaces=table;	%	Phase-to-phase interfaces
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

			[Segmts,OuterLoop,InnerLoops,G.SingularPoints,G.V]=computeSegments(grains);
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