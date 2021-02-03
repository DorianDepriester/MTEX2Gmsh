function [segment_sequences,outLoop,inLoops,sp] = computeSegments(grains)
%COMPUTESEGMENTS Compute the singular point-to-triple point segments defining
% each grain.
%
% [s,out,int,sp,V]=COMPUTESEGMENTS(grains, tol) returns:
%	- s: a list a segments, defined by the nodes' indices
%	- out: the list of outer loops, defined by the segments' indices
%	- int: the list of inner loops, defined by the segments' indices
%	- sp: indices of special points (triple junctions, corners etc.)
% tol is used as tolerance for tracking corners of RoI.
%
%	See also gmshGeo, singularPoints

sp=singularPoints(grains);
sp_all=vertcat(sp{:});	% Store all singular points in a single array

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
        sp_loc=sp_all(ismember(sp_all,poly));
        if ~isempty(sp_loc)
            if poly(end)==poly(1)
                poly(end)=[];                   % 'Open' the loop before shifting
            end
            id_start=find(poly==sp_loc(1),1,'first');
            poly=circshift(poly,[-id_start+1,0]);	% Starts from a special point
        end
        n_tp=length(sp_loc);
        segments=cell(n_tp+1,1);        
        if poly(end)~=poly(1)
            poly(end+1)=poly(1);                %#ok<AGROW> % (re)close the loop
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
vtx_IDs = Vtx_onBounds(grains.boundary);	% IDs of the vertices on the domain boundaries

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

%% Store singular points as structure
sp=cell2struct(vertcat(sp,{[]}),{'triplePoints','quadruplePoints','corners','doublePointsOnBorder','symmetric'});

close(h)

    
function [ segments,idx ] = addSequence(segments,Seq,phaseID)
%Look for existing sequence into a cell array. Add the sequence if it does
%not exists yet. The function returns the cell index of the sequence.
	if isempty(Seq)
        idx=[];
        return
	end
	if Seq(1) == Seq(end)	% Ensure that self loops are unique
		[~,idmin]=min(Seq);
		if idmin~=1
			Seq_ordered=circshift(Seq(1:end-1),[-idmin+1 0]);	% Reorder self loop starting from low vertex index
			Seq=[Seq_ordered; Seq_ordered(1)];
		end
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
	
	

