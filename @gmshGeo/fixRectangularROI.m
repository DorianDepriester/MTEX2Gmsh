function G=fixRectangularROI(obj)
%FIXRECTANGULARROI moves some vertices so that the border of the ROI
%becomes rectangular.
%	FIXRECTANGULARROI(obj) first computes the smallest rectangular box
%	bounding all the geometry, then moves the outer vertices to the closest
%	edge of the box. It also move 4 vertices so that the ROI has sharp
%	corners.
%
%	see also simplify, plotSingularPoints
	G=obj;
	border=G.Interfaces{'ROI Border',1};
	border_segs=G.Segments(border{1});
	ids=unique(vertcat(border_segs{:}));
	V=G.V(ids,:);
	minmax=[min(V); max(V)];

	%% Move vertices to border edges
	d=[V(:,1)-minmax(1) V(:,1)-minmax(2) V(:,2)-minmax(3) V(:,2)-minmax(4)];
	[~,a]=min(abs(d),[],2);
	for k=1:4
		[~,j]=ind2sub([2,2],k);
		V(a==k,j)=minmax(k);
	end

	%% Ensure sharp corners
	corners=zeros(2,2);
	V=obj.V;
	for i=1:2
		for j=1:2
			c=[minmax(i,1) minmax(j,2)];
			d2=(V(:,1)-c(1)).^2+(V(:,2)-c(2)).^2;
			[~,id_dmin]=min(d2);
			corners(i,j)=id_dmin;
			V(id_dmin,:)=c;
		end
	end

	%% Update vertices and singular points
	G.V=V;
	G.SingularPoints.corners=corners(:);
	
	%% Consider double points at corners as corners only
	dpb=G.SingularPoints.doublePointsOnBorder;
	G.SingularPoints.doublePointsOnBorder=dpb(~ismember(dpb,corners(:)));
end
