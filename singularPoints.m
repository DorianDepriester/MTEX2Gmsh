function [sp,V] = singularPoints(grains)
%%SINGULARPOINTS List of singular points in the grains
%
%	SINGULARPOINTS(grains) returns an array of indices of all singular
%	points in the grains. The singular points can be either:
%		-triple junctions,
%		-quadruple junctions (or of even higher order),
%		-corners of the ROI,
%		-double points at the boundaries of the ROI.
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
	minmax=[min(V); max(V)];
	corners=zeros(2,2);
	for i=1:2
		for j=1:2
			c=[minmax(i,1) minmax(j,2)];
			d2=(V(:,1)-c(1)).^2+(V(:,2)-c(2)).^2;
			[~,id_dmin]=min(d2);
			corners(i,j)=id_dmin;
			V(id_dmin,:)=c;
		end
	end
		
	%% Double points on boundary
	% Points shared with larger number of edges than number of grains AND shared between at leat two nodes
	itP=full(sum(I_VG,2)< sum(I_VF,2)) & full(sum(I_VG,2)>1);
	dpb=find(itP);
	
	sp={tp; qp; corners(:); dpb};				% Special points

end

