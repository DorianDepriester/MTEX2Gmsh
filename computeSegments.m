function [segment_sequences,outLoop,inLoops,V,sp] = computeSegments(grains)
%COMPUTESEGMENTS Compute the triple point-to-triple point segments defining
% each grain.
%
% [s,out,in,V,TJ]=COMPUTESEGMENTS(grains) returns:
%	- s: a list a segments, defined by the nodes' indices
%	- out: the list of outer loops, defined by the segments' indices
%	- int: the list of inner loops, defined by the segments' indices
%	- V: the coordinates of the nodes
%	- sp: indices of special points (triple junctions, corners etc.)
%
%	See also gmshGeo

tp=grains.triplePoints.id;              % Triple points
qp=calcQuadruplePoints(grains);         % Quadruple points
bp=doublePtOnBound(grains);             % Double points on the boundaries of the ROI
[V,corners]=sharp_corners(grains);		% Ensure sharp corners

sp=[tp; qp; bp; corners];				% Special points

ng=length(grains);
partitions=cell(ng,1);	% List of sections (ie. grains)
h=waitbar(0,'Computing sections of grain boundaries','Name','gmshGeo class');
n_seg=0;				% Number of segments found

%% Choose the best data type for storing the sequences
maxVtx=max(grains.boundary.F(:));
if maxVtx<intmax('uint8')
	datatype='uint8';
elseif maxVtx<intmax('uint16')
	datatype='uint16';
elseif maxVtx<intmax('uint32')
	datatype='uint32';
else
	datatype='uint64';
end

for ig=1:ng
	waitbar(ig/(2*ng),h);
    F=grains(ig).boundary.F;
    polys=EulerPath(F,datatype);	% Compute outer an inner loops
    nb_loops=size(polys,1);
    surface=cell(nb_loops,1);
    for id_bound=1:nb_loops
        poly=polys{id_bound};
        sp_loc=sp(ismember(sp,poly));
        if ~isempty(sp_loc)
            if poly(end)==poly(1)
                poly(end)=[];                   % 'Open' the loop before shifting
            end
            id_start=find(poly==sp_loc(1),1,'first');
            poly=circshift(poly,-id_start+1);	% Starts from a special point
        end
        n_tp=length(sp_loc);
        segments=cell(n_tp+1,1);        
        if poly(end)~=poly(1)
            poly(end+1)=poly(1);                % (re)close the loop
        end
        Seq=zeros(length(poly),1,datatype);
        id_seq=0;
        id_pt=1;
        id_poly=1;
        if ~isempty(sp_loc)	% poly contains one (or more) special point
            while id_poly<=length(poly)
                Seq(id_pt)=poly(id_poly);
                if id_pt~=1 && any(poly(id_poly)==sp_loc)	% Stops if a special point is reached...
                    Seq((id_pt+1):end)=[];
                    id_seq=id_seq+1;
                    segments{id_seq}=Seq;
					n_seg=n_seg+1;
                    Seq=zeros(length(sp_loc),1,datatype);
                    id_pt=1;
				else										% ...keep on otherwise
                    id_pt=id_pt+1;
                    id_poly=id_poly+1;
                end
            end
		else				% One single sequence
            id_seq=1;
            segments{id_seq}=poly;
			n_seg=n_seg+1;
        end
        surface{id_bound}=segments(1:id_seq);    
    end
    partitions{ig}=surface;
end

%% Indexing all segments
waitbar(0.5,h,'Indexing the segments and looking for duplicates');
segment_sequences=cell(n_seg,2);
vtx_IDs = Vtx_onBounds(grains);	% IDs of the vertices on the domain boundaries

inLoops=cell(ng,1);
outLoop=cell(ng,1);
for i=1:ng
	F=grains(i).boundary.F;
	phaseIds=grains(i).boundary.phaseId;
	waitbar(0.5+i/(2*ng),h);
    parti=partitions{i};
	nloop=length(parti);
	inLoops{i}=cell(nloop-1,1);
	for j=1:nloop
        loopi=parti{j};
        id_segments=zeros(length(loopi),1,'int16');
		for k=1:length(loopi)
			segment=loopi{k};
			phaseID=phaseIds(find( segment(1)==F(:,1) & segment(2)==F(:,2) | segment(2)==F(:,1) & segment(2)==F(:,2) | segment(1)==F(:,end) & segment(2)==F(:,1) | segment(2)==F(:,end) & segment(2)==F(:,1),1,'first'),:);
			if all(ismember(segment,vtx_IDs))	% If all the vertices are on the domain boundary
				segment=segment([1 end]);		% Just keep the end-points
			end
            [segment_sequences,id_segments(k)]=addSequence(segment_sequences,segment,phaseID);
		end
		if j==1
			outLoop{i}=id_segments(:);
		else
			inLoops{i}{j-1}=id_segments(:);
		end
	end
end

% Remove empty segments
lengths=cellfun('length',segment_sequences(:,1));
ide=find(lengths~=0,1,'last');
segment_sequences=segment_sequences(1:ide,:);
close(h)

    
function [ segments,idx ] = addSequence(segments,Seq,phaseID)
%Look for existing sequence into a cell array. Add the sequence if it does
%not exists yet. The function returns the cell index of the sequence.
	if isempty(Seq)
        idx=[];
        return
	end
	lengths=cellfun('length',segments(:,1));
	if all(lengths==0)
		idx=1;
		segments{idx,1}=Seq;
		segments{idx,2}=phaseID;
		return
	end
	idx=0;
	eqlength=find(lengths==length(Seq));
	for i=1:length(eqlength)
		idx_loc=eqlength(i);
		seg_idx=segments{idx_loc,1};
		if isequal(Seq,seg_idx)
			idx=idx_loc;
			return
		elseif isequal(flipud(Seq),seg_idx)
			idx=-idx_loc;
			return
		elseif Seq(1)==Seq(end) && seg_idx(1)==seg_idx(end) && all(unique(Seq)==unique(seg_idx))
			idx=idx_loc;
			return
		end
	end
	if idx==0
		idx=find(lengths==0,1,'first');
		if isempty(idx)
			idx=length(segments)+1;
		end
		segments{idx,1}=Seq;
		segments{idx,2}=phaseID;
	end
	
function qp = calcQuadruplePoints(grains)
    gB=grains.boundary;
    I_VF = gB.I_VF;
    I_VG = (I_VF * gB.I_FG)==2;
    itP = full(sum(I_VG,2)>2 & sum(I_VF,2)>3);	% Due to Voronoi decomposition, the vertex order can actually be higher than 4
    qp=find(itP);

function bp = doublePtOnBound(grains)
    gB=grains.boundary;
    V=grains.boundary.V;
    I_VF = gB.I_VF;
    I_VG = (I_VF * gB.I_FG)==2;

    Vtx_onBounds_Ids=Vtx_onBounds(grains);      % Corresponding vertices
    on_bound=false(length(V),1);                % Vertex on bound or not?
    on_bound(Vtx_onBounds_Ids(:))=1;
    dp=full(sum(I_VG,2)>1 & sum(I_VF,2)>1);     % Double points
    dp_onBound = dp & on_bound;
    
    bp=find(dp_onBound);
	
function [ vtx_IDs ] = Vtx_onBounds(grains)
    gB=grains.boundary;
    outerBoundary_id=any(gB.grainId==0,2);      % IDs of the faces neighbouring no other grains
    list_IDs=gB.F(outerBoundary_id(:),:);       % Corresponding vertices
    vtx_IDs=unique(list_IDs);

function [V,corners]=sharp_corners(grains)
    gB=grains.boundary;
    outerBoundary_id=any(gB.grainId==0,2);      % IDs of the faces neighbouring no other grains
    list_IDs=gB.F(outerBoundary_id(:),:);       % Corresponding vertices
    vtx_IDs=unique(list_IDs);
	V=gB.V;
    Vl=V(vtx_IDs,:);                          % Coordinates of the vertices on outer boundary

    x_lim= [min(gB.V(:,1)) max(gB.V(:,1))];
    y_lim= [min(gB.V(:,2)) max(gB.V(:,2))]; 
    lim= [x_lim y_lim];

    %% Move outer vertices toward the closest limits
    for i=1:length(vtx_IDs)
        dist=[(Vl(i,1)-lim(1)).^2 (Vl(i,1)-lim(2)).^2 (Vl(i,2)-lim(3)).^2 (Vl(i,2)-lim(4)).^2];
        [~,I]=min(dist);
        if I==1 || I==2
            Vl(i,1)=lim(I);
        else
            Vl(i,2)=lim(I);
        end
    end

    %% Ensure that an end-point is located in each corner
    corners=zeros(2,2); % IDs of the vertices forming the 4 corners
    for i=1:2
        for j=1:2
            dist=(Vl(:,1)-x_lim(i)).^2+(Vl(:,2)-y_lim(j)).^2;
            [~,I]=min(dist);
            Vl(I,:)=[x_lim(i),y_lim(j)];
            corners(i,j)=vtx_IDs(I);
        end
    end
    corners=corners(:);
    
    V(vtx_IDs,:)=Vl; % Update the coordinates of the vertices
