function segList = borderLoop(G)
	Border=G.Interfaces{'ROI Border',1};
	Border=Border{1};
	F=vertcat(G.Segments{Border});
	F=[F(1:2:end) F(2:2:end)];
	dataType='single';
	p=EulerPath(F,dataType);
	p=p{1};
	nseg=length(p)-1;
	segList=zeros(nseg,1,dataType);			%	Cast the segment list like that of other loops
	for i=1:nseg
		I=find(F(:,1)==p(i) & F(:,2)==p(i+1),1,'first');
		if isempty(I)
			I=find(F(:,2)==p(i) & F(:,1)==p(i+1),1,'first');
			segList(i)=-cast(Border(I),dataType);	%	Border(I) is unsigned 
		else
			segList(i)=Border(I);
		end
	end
	segList=cast(segList,'like',G.Grains.OuterLoop{1});
end

