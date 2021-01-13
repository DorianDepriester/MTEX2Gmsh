function s=evalElementSize(obj)
	%EVALELEMENTSIZE Automatically evaluate the element size.
	%	EVALELEMENTSIZE(Object) computes the mean node-to-node
	%	distance in each segment.
	%
	%	Note: the value from EVALELEMENTSIZE is used by default in
	%	the savegeo method.
	%
	%	See also savegeo.
	segmts=obj.Segments;
	nseg=length(segmts);
	d=zeros(nseg,1);
	wt=zeros(nseg,1);
	for i=1:nseg
		if length(segmts{i})>2
			XY=obj.V(segmts{i},:);
			dl=sum((XY(2:end,:)-XY(1:end-1,:)).^2,2);
			d(i)=mean(sqrt(dl));
			wt(i)=length(dl);
		end
	end
	keep= d~=0;
	d=d(keep);
	wt=wt(keep);
	s=d'*wt/sum(wt);
end