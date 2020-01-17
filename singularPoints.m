function sp = singularPoints(grains)
%%SINGULARPOINTS List of singular points in the grains
%
%	SINGULARPOINTS(grains) returns an array of indices of all singular
%	points in the grains. The singular points can be either:
%		-triple junctions,
%		-quadruple junctions (or of even higher order),
%		-points at the boundaries of the ROI.
%
%	See also cond_smooth

	%% Triple points
	tp=grains.triplePoints.id;
	
	%% Quadruple points
	gB=grains.boundary;
    I_VF = gB.I_VF;
    I_VG = (I_VF * gB.I_FG)==2;
    itP = full(sum(I_VG,2)>2 & sum(I_VF,2)>3);	% Due to Voronoi decomposition, the vertex order can actually be higher than 4
    qp=find(itP);
	
	%% Points on the boundary of the ROI
    V=grains.boundary.V;
    gB=grains.boundary;
    outerBoundary_id=any(gB.grainId==0,2);      % IDs of the faces neighbouring no other grains
    list_IDs=gB.F(outerBoundary_id(:),:);       % Corresponding vertices
    Vtx_onBounds_Ids=unique(list_IDs);
    on_bound=false(length(V),1);                % Vertex on bound or not?
    on_bound(Vtx_onBounds_Ids(:))=1;
    bp=find(on_bound);
	
	sp=[tp; qp; bp];				% Special points

end

