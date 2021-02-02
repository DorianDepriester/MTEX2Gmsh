function [Gr,varargout]=addSymmetry(obj,direction)
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
%	[G,A_seg]=ADDSYMMETRY(G) also returns the adjancency matrix for
%	segments (giving indices segments symmetrical to each other).
%
%	[G,A_seg,A_V]=ADDSYMMETRY(G) does the same, plus it gives
%	correspondance between vertices the same way.
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
	n_oldV=size(V,1);
	
	%% Track if vertex has be symmetrized
	Vsym=NaN(size(Vdpb));			% List of new vertices
	A_V=zeros(size(Vdpb,1),2);	% Adjancency table
	n_newV=0;
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
			n_newV=n_newV+1;
			A_V(n_newV,1)=dpb(i);
			d2=(Vdpb(:,1)-new_coords(1)).^2+(Vdpb(:,2)-new_coords(2)).^2;
			duplicate= d2<=tol^2;
			if all(~duplicate)
				% No duplicate found
				Vsym(n_newV,:)=new_coords;
				A_V(n_newV,2)=n_newV+n_oldV;
			else
				% Slightly change the coordinate of the existing vertex so 
				% that it is exactly symmetric
				Vdpb(duplicate,1)=new_coords(1);
				Vdpb(duplicate,2)=new_coords(2);
				A_V(n_newV,2)=dpb(duplicate);
			end
		end
	end
	Vsym=Vsym(1:n_newV,:);		% New vertices due to symmetry
	A_V=A_V(1:n_newV,:);
	

	%% Split segments where vertices are added
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
			% Check if the new vertex is in between two existing vertices
			if abs(det(M))<1e-10 && (proj < dot(v,v)) && (proj > 0)
				id_V=i+n_oldV;
				Segments{id_seg}=[seg_j(1); id_V];	% Split existing segment
				id_new_seg=length(Segments)+1;
				Segments{id_new_seg}=[id_V; seg_j(2)];	% Add new segment
				border=[border; id_new_seg];			%#ok<AGROW> % Add reference to new segment to the border list
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
	
	%% Compute adjacency between segments (if needed)
	if nargout>1
		n_seg=length(border);
		A_seg=zeros(n_seg,2);
		n_seg_adj=0;
		for i=1:n_seg
			border_i=border(i);
			ref_Vs=Segments{border_i};
			oppo_seg=zeros(2,1);
			for k=1:2
				j=find(A_V==ref_Vs(k));
				if isempty(j)
					break
				end
				[I,J]=ind2sub(size(A_V),j);
				oppo_seg(k)=A_V(I,3-J);
			end
			if all(oppo_seg)
				if ~ismember(border_i,A_seg);
					n_seg_adj=n_seg_adj+1;
					A_seg(n_seg_adj,1)=border_i;
					a= cellfun(@(x) isequal(x,oppo_seg) | isequal(x,flip(oppo_seg)),Segments);
					A_seg(n_seg_adj,2)=find(a);
				end
			end
		end
		varargout{1}=A_seg(1:n_seg_adj,:);
	end
	
	%% Update properties for the output object
	Gr.Interfaces{'ROI Border',1}={border};
	Gr.V=V;
	Gr.V(dpb,:)=Vdpb;
	Gr.Grains.OuterLoop=outerLoops;
	Gr.Segments=Segments;
	id_sym=(1:n_newV)'+n_oldV;	% Ids for all the new vertices
	Gr.SingularPoints.Symmetric=cast(id_sym,'like',Segments{1});	% Create a new property for singular points
	
	if nargout==3
		varargout{2}=A_V;
	end
end
	