function obj = regularize(obj, varargin)
%REGULARIZE Collapse the segments shorter than a given threshold into vertices.
%	OBJ = REGULARIZE(OBJ) removes every segment whose curvilinear length is 
%   lower than a threshold value. The endpoints of each removed segment are
%	merged into a single vertex located at the barycenter of its knots; all
%	the remaining segments, the grain loops, the interfaces and the
%	singular points are updated accordingly. The default threshold value is
%	set to the spatial resolution.
%
%	Collapsible segments that are adjacent (sharing an endpoint) are merged
%	into a single vertex, so that their shared singular point does not
%	split. This also holds at quadruple points, where more than two
%	neighboring segments meet.
%
%	Segments lying on the border of the ROI are never collapsed, and a
%	short segment with one endpoint on the ROI border is collapsed onto
%	that endpoint, so that the border of the domain (and its corners)
%	remains strictly fixed. Loops that would be left with less than three
%	segments are left unchanged.
%
%	Note that the segment IDs stored in the loops are signed: a negative
%	ID means that the segment is traversed against its stored knot
%	sequence.
%
%   REGULARIZE(G) use the spatial resolution of gmshGeo object G as a
%   threshold value.
%
%   REGULARIZE(G, epsilon) use EPSILON as a threshold value.
%
%	See also gmshGeo, gmshGeo/mesh, fixRectangularROI.

    ip = inputParser;
    addOptional(ip,'epsilon', obj.Resolution)
    parse(ip,varargin{:});
    epsilon = ip.Results.epsilon;

	seg = obj.Segments;
	n_seg = length(seg);
	L = obj.segmentLenghts;
	short_segment = find(L < epsilon).';
	if isempty(short_segment)
		return
	end

	grains = obj.Grains;
	outerLoops = grains.OuterLoop;
	innerLoops = grains.InnerLoops;

	% Segments lying on the ROI border must be preserved
	borderSegs = [];
	row = strcmp(obj.Interfaces.Properties.RowNames, 'ROI Border');
	if any(row)
		borderSegs = double(obj.Interfaces.SegmentIDs{row}(:)).';
	end
	borderVtx = [];
	if ~isempty(borderSegs)
        a = cellfun(@(s) double(s(:)), seg(borderSegs), 'UniformOutput', false);
		borderVtx = unique(vertcat(a{:}));
	end
	collapsible = setdiff(short_segment, borderSegs);

	% Protect the loops that would be left with less than 3 segments
	for i = 1:height(grains)
		loops = [{outerLoops{i}}; innerLoops{i}];
		for j = 1:length(loops)
			loop = abs(double(loops{j}(:)));
			n_short = sum(ismember(loop, collapsible));
			if length(loop) - n_short < 3
				collapsible = setdiff(collapsible, loop);
			end
		end
	end

	% Protect the segments with both endpoints on the ROI border
	for k = collapsible
		s = double(seg{k}(:));
		if all(ismember(s([1 end]), borderVtx))
			collapsible = setdiff(collapsible, k);
		end
	end
	if isempty(collapsible)
		return
	end

	% Group the collapsible segments into chains of adjacent segments
	% (sharing an endpoint): all the segments of a chain are merged into
	% a single vertex
	n_col = length(collapsible);
	parent = (1:n_col).';
	for k = 1:n_col
		s = double(seg{collapsible(k)}(:));
		for v = s([1 end]).'
			for m = (k+1):n_col
				sm = double(seg{collapsible(m)}(:));
				if any(sm([1 end]) == v)
					parent = unionSeg(parent, k, m);
				end
			end
		end
	end
	root = zeros(n_col, 1);
	for k = 1:n_col
		root(k) = findRoot(parent, k);
	end

	% vmap maps each endpoint of a collapsed segment to the vertex of its
	% chain: the barycenter of the chain knots, or the chain endpoint
	% lying on the ROI border if any
    V = obj.V;
	n_vtx = size(V, 1);
	vmap = zeros(n_vtx, 1);
	new_vertices = zeros(0, size(V, 2));
	for r = unique(root).'
		members = find(root == r);
		ends = [];
		for m = members.'
			sm = double(seg{collapsible(m)}(:));
			ends = [ends; sm([1 end])]; %#ok<AGROW>
		end
		ends = unique(ends);
		b = ends(ismember(ends, borderVtx));
		if ~isempty(b)
			mid = b(1);
        else
            a = cellfun(@(m) double(seg{collapsible(m)}(:)), num2cell(members), 'UniformOutput', false);
			knots = unique(vertcat(a{:}));
			new_vertices = [new_vertices; mean(V(knots, :), 1)]; %#ok<AGROW>
			mid = n_vtx + size(new_vertices, 1);
		end
		vmap(ends) = mid;
	end

	% Remove the collapsible segments from the loops
	for i = 1:height(grains)
		outerLoops{i} = outerLoops{i}(~ismember(abs(outerLoops{i}), collapsible));
		inner = innerLoops{i};
		for j = 1:length(inner)
			inner{j} = inner{j}(~ismember(abs(inner{j}), collapsible));
		end
		innerLoops{i} = inner;
	end

	% Replace the collapsed endpoints in the endpoints of all the segments
	deg = false(n_seg, 1);
	for j = 1:n_seg
		k = double(seg{j}(:));
		deg(j) = vmap(k(1)) ~= 0 && vmap(k(end)) ~= 0 && vmap(k(1)) == vmap(k(end));
		if vmap(k(1)) ~= 0
			k(1) = vmap(k(1));
		end
		if vmap(k(end)) ~= 0
			k(end) = vmap(k(end));
		end
		seg{j} = k;
	end

	% Remove the collapsed (and degenerate) segments and remap all the IDs
	keep = true(n_seg, 1);
	keep(collapsible) = false;
	keep(deg) = false;
	remap = zeros(n_seg, 1);
	remap(keep) = 1:sum(keep);
	seg = seg(keep);
	for i = 1:height(grains)
		loop = sign(double(outerLoops{i}(:))) .* remap(abs(double(outerLoops{i}(:))));
		outerLoops{i} = loop(loop ~= 0);
		inner = innerLoops{i};
		for j = 1:length(inner)
			inner{j} = sign(double(inner{j}(:))) .* remap(abs(double(inner{j}(:))));
			inner{j}(inner{j} == 0) = [];
		end
		innerLoops{i} = inner;
	end
	interfaces = obj.Interfaces;
	for i = 1:height(interfaces)
		ids = remap(interfaces.SegmentIDs{i});
		interfaces.SegmentIDs{i} = ids(ids > 0);
	end
	sp = obj.SingularPoints;
	field_names = fieldnames(sp);
	for i = 1:length(field_names)
		ids = double(sp.(field_names{i})(:));
		m = vmap(ids);
		ids(m ~= 0) = m(m ~= 0);
		sp.(field_names{i}) = ids;
	end

	obj.Segments = seg;
	obj.V = [V; new_vertices];
	grains.OuterLoop = outerLoops;
	grains.InnerLoops = innerLoops;
	obj.Grains = grains;
	obj.Interfaces = interfaces;
	obj.SingularPoints = sp;
end

function p = unionSeg(p, a, b)
	ra = findRoot(p, a);
	rb = findRoot(p, b);
	if ra ~= rb
		p(rb) = ra;
	end
end

function r = findRoot(p, x)
	r = x;
	while p(r) ~= r
		r = p(r);
	end
end