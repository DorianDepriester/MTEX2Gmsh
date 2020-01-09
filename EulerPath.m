function polys=EulerPath(F,datatype)
%EULERPATH Compute the Euler path of a graph.
%
%	EULERPATH(F,DTYPE) where F is an N-2 array defining the connectivity
%	between nodes (defined by their indices). The return value is a cell 
%	array defined by nodes indices written in the DTYPE datatype (can be
%	'uint8', 'uint16', 'uint32' or 'uint64').
%
%	See also computeSegments

vtx_list=unique(F);
mult=hist(single(F(:)),single(vtx_list))';
nbif=nnz(mult>2);	% Number of triple nodes (ambiguous path)
npaths=3^nbif;
innerLoops=cell(npaths,1);
outerLoops=cell(1);
configs=combvec3(nbif);
nf=size(F,1);
visited=false(size(vtx_list));

Loop=zeros(nf,1,datatype);
Inner=false(1,npaths);

connect_id=1;	% id of the connected region (if holes)
while any(~visited)
	start=vtx_list(find(~visited & mult==2,1,'first'));
	Loop(1)=start;
	ids=find(F==start,1,'first');
	[is,js]=ind2sub(size(F),ids);
	if js==1
		Loop(2)=F(is,2);
	else
		Loop(2)=F(is,1);
	end
	for k=1:npaths
		bifurcId=1;		% If the Euler tour is not unique (multiple point)
		vt1=Loop(1);
		vt2=Loop(2);
		visited(vtx_list==vt1)=1;
		visited(vtx_list==vt2)=1;
		config=configs(k,:);
		for i=3:(nf+1)
			j=find(F(:,1)==vt2 & F(:,2)~=vt1 | F(:,2)==vt2 & F(:,1)~=vt1);	% avoid backward walk
			nj=length(j);
			if isempty(j)
				break
			elseif nj>1
				Fj=F(j(config(bifurcId)),:);
				bifurcId=bifurcId+1;
			else 
				Fj=F(j,:);
			end
			Loop(i)=Fj(Fj~=vt2);
			vt1=vt2;
			vt2=Loop(i);
			visited(vtx_list==vt2)=1;
			ip=find(vt2==Loop(1:i-1),1,'first');
			if ~isempty(ip) 
				if ip==1
					outerLoops{connect_id}=Loop(1:i);
				else
					innerLoops{k,connect_id}=Loop(ip:i);
					Inner(k)=1;
				end
				break
			end
		end
	end
	connect_id=connect_id+1;
end

%% Sort loops such that the (real) outer loop appears first
innerLoops=innerLoops(:);
outerLoops=outerLoops(:);
loop_len=cellfun('length',outerLoops);
[~,order]=sort(loop_len,'descend');	% The largest one is considered as the outer loop
outerLoops=outerLoops(order);		% Rearrange with decreasing length
innerLoops=[outerLoops(2:end) ; innerLoops];

%% Remove empty sets and duplicates
for i=1:length(innerLoops)
	Loopi=innerLoops{i};
	for j=1:(i-1)
		if length(Loopi)==length(innerLoops{j}) && ( all(flipud(Loopi)==innerLoops{j}) || all(Loopi==innerLoops{j}) )
			innerLoops{i}=[];
			break
		end
	end
end
% Remove empty loops
innerLoops(cellfun('isempty',innerLoops))=[];
outerLoops(cellfun('isempty',outerLoops))=[];

polys=[outerLoops(1); innerLoops];



function M = combvec3(n)
	M=zeros(3^n,n,'uint32');
	for i=1:n
		M(:,i)=repmat([ones(3^(i-1),1) ; 2*ones(3^(i-1),1) ; 3*ones(3^(i-1),1)],3^(n-i),1);
	end
