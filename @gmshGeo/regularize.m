function obj = regularize(obj, epsilon)
%REGULARIZE Collapse the segments shorter than a given threshold into vertices.
%	OBJ = REGULARIZE(OBJ, EPSILON) removes every segment whose curvilinear
%	length is lower than EPSILON. Each removed segment is replaced by a
%	single vertex located at the barycenter of its knots; the neighbouring
%	segments, the grain loops and the interfaces are updated accordingly.
%	Collapsible segments that are adjacent (sharing an endpoint) are merged
%	into a single vertex, so that their shared triple point does not split.
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
%	See also gmshGeo, mesh, fixRectangularROI, simplify.

	seg = obj.Segments;
	V = obj.V;
	n_seg = length(seg);
	L = zeros(n_seg, 1);
	for i = 1:n_seg
		L(i) = BS_length(V(seg{i}, :));
	end
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
		borderVtx = unique(vertcat(cellfun(@(s) double(s(:)), seg(borderSegs), 'UniformOutput', false){:}));
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

	% Each chain is merged into a single vertex: the barycenter of its
	% knots, or its endpoint if it lies on the ROI border
	merge_ids = zeros(n_col, 1);
	new_vertices = zeros(0, size(V, 2));
	for r = unique(root).'
		members = find(root == r);
		ends = [];
		for m = members.'
			sm = double(seg{collapsible(m)}(:));
			ends = [ends; sm([1 end])]; %#ok<AGROW>
		end
		b = unique(ends(ismember(ends, borderVtx)));
		if ~isempty(b)
			merge_ids(members) = b(1);
		else
			knots = unique(vertcat(cellfun(@(m) double(seg{collapsible(m)}(:)), num2cell(members), 'UniformOutput', false){:}));
			new_vertices = [new_vertices; mean(V(knots, :), 1)]; %#ok<AGROW>
			merge_ids(members) = size(V, 1) + size(new_vertices, 1);
		end
	end

	for i = 1:height(grains)
		[outerLoops{i}, seg] = regularizeLoop(outerLoops{i}, collapsible, merge_ids, seg);
		inner = innerLoops{i};
		for j = 1:length(inner)
			[inner{j}, seg] = regularizeLoop(inner{j}, collapsible, merge_ids, seg);
		end
		innerLoops{i} = inner;
	end

	% Remove the collapsed segments and remap all the IDs
	keep = true(n_seg, 1);
	keep(collapsible) = false;
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

	obj.Segments = seg;
	obj.V = [V; new_vertices];
	grains.OuterLoop = outerLoops;
	grains.InnerLoops = innerLoops;
	obj.Grains = grains;
	obj.Interfaces = interfaces;
end

function [loop, seg] = regularizeLoop(loop, collapsible, merge_ids, seg)
%REGULARIZELOOP Replace each occurrence of the collapsible segments in LOOP
%	by a single vertex spliced into the neighbouring segments.

	k = 1;
	while k <= length(collapsible)
		a = find(abs(loop) == collapsible(k), 1, 'first');
		if isempty(a) || length(loop) <= 3
			k = k + 1;
			continue
		end
		if a == 1
			prev = length(loop);
			next = 2;
		elseif a == length(loop)
			prev = a - 1;
			next = 1;
		else
			prev = a - 1;
			next = a + 1;
		end
		id_new = merge_ids(k);
		id_prev = loop(prev);
		s = seg{abs(id_prev)}(:);
		if id_prev > 0		% Loop end = last stored knot
			seg{abs(id_prev)} = [s(1:end-1); id_new];
		else				% Loop end = first stored knot
			seg{abs(id_prev)} = [id_new; s(2:end)];
		end
		id_next = loop(next);
		s = seg{abs(id_next)}(:);
		if id_next > 0		% Loop start = first stored knot
			seg{abs(id_next)} = [id_new; s(2:end)];
		else				% Loop start = last stored knot
			seg{abs(id_next)} = [s(1:end-1); id_new];
		end
		loop(a) = [];
	end
end

function L = BS_length(V)
	bs = BSpline(V);
	L = sum(vecnorm(diff(bs, 1, 1), 2, 2));
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
