function sp = singularPoints(grains)
%%SINGULARPOINTS List of singular points in the grains
%
%	SINGULARPOINTS(grains,tol) returns an array of indices of all singular
%	points in the grains. The singular points can be either:
%		-triple junctions,
%		-quadruple junctions (or of even higher order),
%		-corners of the ROI,
%		-double points at the boundaries of the ROI.
%	The tolerance (tol) is used to determine if a vertex can be used as a
%	corner for ROI or not.
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
	
	%% Corners of the ROI
	% Find closest points to the corners of ROI, and move them to the
	% corners.
	V=grains.V;
	corners=zeros(2,2);
	minmax=[min(V); max(V)];
	meanV2V=mean(grains.perimeter./cellfun('length',grains.poly));	% Mean vertex-to-vertex distance
	tol=meanV2V/10;
	for i=1:2
		for j=1:2
			c=[minmax(i,1) minmax(j,2)];
			d2=(V(:,1)-c(1)).^2+(V(:,2)-c(2)).^2;
			[d2min,id_dmin]=min(d2);
			if sqrt(d2min)<tol
				corners(i,j)=id_dmin;
				V(id_dmin,:)=c;
			end
		end
	end
	corners=corners(corners~=0);
		
	%% Double points on boundary
	border=Vtx_onBounds(gB);
	% Points shared with more than 2 faces
	itP=full(sum(I_VF,2)>2);
	dpb=intersect(find(itP),border);
	dpb=dpb(~ismember(dpb,corners));	% Avoid double points at the corners
	
	sp={tp; qp; corners; dpb};			% Special points

end

