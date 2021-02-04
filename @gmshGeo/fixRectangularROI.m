function Gr=fixRectangularROI(obj)
%FIXRECTANGULARROI moves some vertices so that the border of the ROI
%becomes rectangular.
%	FIXRECTANGULARROI(obj) first computes the smallest rectangular box
%	bounding all the geometry, then moves the outer vertices to the closest
%	edge of the box. It also moves (or creates) 4 vertices so that the ROI
%	has sharp corners.
%
%	see also simplify, plotSingularPoints
	Gr=obj;
	V=Gr.V;

	tol=Gr.Resolution/10;
	border=Gr.Interfaces{'ROI Border',1};
	border=border{1};
	border_segs=Gr.Segments(border);
	Segments=Gr.Segments;
	minmax=[min(V); max(V)];
	outerLoops=Gr.Grains{:,'OuterLoop'};

	ids_Vb=unique(vertcat(border_segs{:}));



	%% Move vertices to border edges
	d=[V(ids_Vb,1)-minmax(1) V(ids_Vb,1)-minmax(2) V(ids_Vb,2)-minmax(3) V(ids_Vb,2)-minmax(4)];
	[~,a]=min(abs(d),[],2);
	for k=1:4
		[~,j]=ind2sub([2,2],k);
		V(ids_Vb(a==k),j)=minmax(k);
	end

	%% Find/create sharp corners
	corners=zeros(2,2);
	for i=1:2
		for j=1:2
			C=[minmax(i,1) minmax(j,2)];
			% First, try with existing vertices
			d2=(V(ids_Vb,1)-C(1)).^2+(V(ids_Vb,2)-C(2)).^2;
			[d2_min,id_dmin]=min(d2);
			if d2_min<tol^2
				corners(i,j)=ids_Vb(id_dmin);
				V(ids_Vb(id_dmin),:)=C;
			else
				% No candidate vertex found -> need to create a new one
				id_border_min=0;
				dmin=inf;			
				for k=1:length(border_segs)
					% Fetch the best segment to put new vertex in
					seg_k=border_segs{k};
					V1=V(seg_k(1),:);
					V2=V(seg_k(2),:);
					u=V2-V1;
					v=C-V1;
					u2=sum(u.^2);
					d=abs(det([u; v]))/sqrt(u2);
					proj=dot(u,v);
					if d<dmin && (proj < u2) && (proj > 0)
						id_border_min=k;
						dmin=d;
					end
				end
				V=[V; C];	% Add new vertex
				id_new_V=size(V,1);
				corners(i,j)=id_new_V;
				id_seg_min=border(id_border_min);
				old_seg=Segments{id_seg_min};
				Segments{id_seg_min}=[old_seg(1); id_new_V];	% Split the segment
				new_seg=[id_new_V; old_seg(2)];					% Create a new one
				Segments{end+1}=new_seg;
				id_new_seg=length(Segments);
				for k=1:length(outerLoops)
					if ismember(id_seg_min,outerLoops{k})
						outerLoops{k}=[outerLoops{k}; id_new_seg];	% Update definition for the grains
					end
				end
				border=[border; id_new_seg];	% Add the new segment to border list		
			end
		end
	end

	Gr.Interfaces{'ROI Border',1}={border};
	Gr.V=V;
	Gr.Grains.OuterLoop=outerLoops;
	Gr.Segments=Segments;
	Gr.SingularPoints.corners=corners(:);
end
