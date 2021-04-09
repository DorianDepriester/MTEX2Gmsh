function [LineLoops,PlaneSurface]=uniqueLoops(Grains)
	h=waitbar(0,'Removing loop duplicates','Name','Closed loops');
	n_grains=height(Grains);
	PlaneSurface=cell(height(Grains),1);
	nLoops=n_grains+sum(cellfun(@length,Grains.InnerLoops));	%	Overall number of loops
	LineLoops=cell(nLoops,1);
	jmax=0;	%	Number of unique loops
	for i=1:n_grains
		waitbar(i/n_grains,h);
		OuterLoop=Grains.OuterLoop{i};
		InnerLoops=Grains.InnerLoops{i};
		Seq=zeros(1+length(InnerLoops),1,'int32');
		new=true;
		for j=1:size(LineLoops,1)
			if isempty(LineLoops{j})
				break
			end
			if isequal(sort(abs(OuterLoop)),sort(abs(LineLoops{j})))
				new=false;
				break
			end
		end
		if new
			LineLoops{j}=OuterLoop;
		end
		Seq(1)=j;
		jmax=max(j,jmax);

		for k=1:length(InnerLoops)
			new=true;			
			for j=1:size(LineLoops,1)
				if isempty(LineLoops{j})
					break
				end
				if isequal(sort(abs(InnerLoops{k})),sort(abs(LineLoops{j})))
					new=false;
					break
				end
			end
			if new
				LineLoops{j}=InnerLoops{k};
			end
			Seq(1+k)=j;
			jmax=max(j,jmax);			
		end
		PlaneSurface{i}=Seq;
	end
	LineLoops=LineLoops(1:jmax);
	close(h);
end

