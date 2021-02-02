function Gr=addSymmetry(obj,direction)
%ADDSYMMETRY adds symmetry conditions on the geometry, so that singular
%points exist are on opposite sides of the RoI. This function is intended
%to mesh the geometry with periodic conditions.
%
%	ADDSYMMETRY(G,'x') inserts dummy singular points at the borders of the
%	ROI so that those borders are symmetric along the X axis (i.e.
%	reflection with respect to the median axis parallel to Y).
%
%	ADDSYMMETRY(G,'y') does the same along the Y axis.
%
%	ADDSYMMETRY(G,'both') applies symmetry along both X and Y axes.
%
%	NOTE: this function only works on perfectly rectangular RoI. Thus, the
%	fixRectangularROI function is applied on the geometry.
%
%	see also fixRectangularROI

	xsym=strcmpi(direction,'x') || strcmpi(direction,'both');
	ysym=strcmpi(direction,'y') || strcmpi(direction,'both');
		
	Gr=fixRectangularROI(obj);
	dpb=Gr.SingularPoints.doublePointsOnBorder;
	V=Gr.V;
	Vdpb=V(dpb,:);			% Singular point to be symmetrized

	npt=length(dpb);
	minmax=[min(Vdpb); max(Vdpb)];
	tol=Gr.Resolution/10;
	
	%% Track if vertex has be symmetrized
	Vsym=NaN(size(Vdpb));	% List of new vertices
	for i=1:npt
		new_coords=[];
		if xsym
			if abs(Vdpb(i,1)-minmax(1,1))<tol			% Vertex close to left border
				new_coords=[minmax(2,1) Vdpb(i,2)];
			elseif abs(Vdpb(i,1)-minmax(2,1))<tol		% Vertex close to right border
				new_coords=[minmax(1,1) Vdpb(i,2)];
			end
		end
		if ysym
			if abs(Vdpb(i,2)-minmax(1,2))<tol			% Vertex close to bottom border
				new_coords=[Vdpb(i,1) minmax(2,2)];
			elseif abs(Vdpb(i,2)-minmax(2,2))<tol
				new_coords=[Vdpb(i,1) minmax(1,2)];		% Vertex close to top border
			end
		end
		if ~isempty(new_coords)
			% Avoid duplicates: check that the vertices does not exist yet
			d2=(Vdpb(:,1)-new_coords(1)).^2+(Vdpb(:,2)-new_coords(2)).^2;
			duplicate= d2<=tol^2;
			if all(~duplicate)
				% No duplicate found
				Vsym(i,:)=new_coords;
			else
				% Slightly change the coordinate of the existing vertex so 
				% that it is exactly symmetric
				Vdpb(duplicate,1)=new_coords(1);
				Vdpb(duplicate,2)=new_coords(2);
			end
		end
	end
	Vsym=Vsym(~isnan(Vsym(:,1)),:);	% New vertices due to symmetry

	%% Split segments where vertices are added
	n_newV=size(Vsym,1);
	n_oldV=size(V,1);
	V=[V; Vsym];
	border=Gr.Interfaces{'ROI Border',1};
	border=border{1};
	Segments=Gr.Segments;
	outerLoops=Gr.Grains{:,'OuterLoop'};
	for i=1:n_newV
		Vi=Vsym(i,:);
		j=1;
		while j<=length(border)
			id_seg=border(j);
			seg_j=Segments{id_seg};
			A=V(seg_j(1),:);
			u=Vi-A;
			v=V(seg_j(2),:)-A;
			M=[u; v];
			proj=dot(u,v);
			% Check if the nex vertex is in between two existing vertices
			if abs(det(M))<1e-10 && (proj < dot(v,v)) && (proj > 0)
				id_V=i+n_oldV;
				Segments{id_seg}=[seg_j(1); id_V];	% Split existing segment
				id_new_seg=length(Segments)+1;
				Segments{id_new_seg}=[id_V; seg_j(2)];	% Add new segment
				border=[border; id_new_seg];			% Add reference to new segment to the border list
				for k=1:length(outerLoops)
					if ismember(id_seg,outerLoops{k})
						outerLoops{k}=[outerLoops{k}; id_new_seg];
					end
				end
				break
			end
			j=j+1;
		end		
	end
	
	%% Update properties for the output object
	Gr.Interfaces{'ROI Border',1}={border};
	Gr.V=V;
	Gr.V(dpb,:)=Vdpb;
	Gr.Grains.OuterLoop=outerLoops;
	Gr.Segments=Segments;
	id_sym=(1:n_newV)'+n_oldV;	% Ids for all the new vertices
	Gr.SingularPoints.Symmetric=cast(id_sym,'like',Segments{1});	% Create a new property for singular points
end
	