function G=simplify(obj,varargin)
	%SIMPLIFY Applies the Douglas-Peucker algorithm to reduce the number of
	%points in each segment.
	%
	%	SIMPLIFY(Object) reduces the number elements with penalty 
	%	lenght equal to one tenth of the default element size.
	%
	%	SIMPLIFY(Object,epsilon) uses epsilon as the penalty
	%	length.
	%
	%	See also plot.
	G=obj;
	if nargin==1
		epsilon=obj.Resolution/10;
	else
		epsilon=varargin{1};
	end

	%% Apply the Douglas-Peucker algorithm on each segment
	segments=obj.Segments;
	h=waitbar(0,'Applying the Douglas-Peucker algorithm','Name','Simplification of the boundaries','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
	setappdata(h,'canceling',0)
	n_seg=length(segments);
	for i=1:n_seg
		if getappdata(h,'canceling')
			delete(h)
			return
		end
		waitbar(i/n_seg,h);
		segmt=obj.Segments{i};
		remains=DouglasPeucker(obj.V(segmt,:),epsilon);
		segments{i}=segmt(remains);
	end

	%% Remove unused vertices and update the segments
	Vtx=G.V;
	used=ismember(1:size(Vtx,1),cell2mat(segments));	% Track whether the vertices are used or not in the segments
	new_idx=cast(cumsum(used),'like',segments{1});
	for j=1:length(segments)
		segments{j}=new_idx(segments{j})';
	end
	G.Segments=segments;			% Update the segments
	G.V=Vtx(used,:);				% Remove unused vertices
	
	%% Update singular points
	sp=G.SingularPoints;
	field_names=fieldnames(sp);
	for i=1:length(field_names)
		sp_i_old=sp.(field_names{i});
		G.SingularPoints.(field_names{i})=new_idx(sp_i_old);	% Update singular points
	end

	
	delete(h)
end
		
function remains=DouglasPeucker(V,epsilon)
	npt=size(V,1);
	if npt>2
		X=V(:,1);
		Y=V(:,2);		
		Pdg=zeros(npt,1);
		Pdg(1)=1;
		Pdg(end)=npt;
		if X(1) == X(end) && Y(1)==Y(end)	%	Closed loop
			mid=round(npt/2);				%	Add mid-point
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