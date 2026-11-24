function obj = regularize(obj, epsilon)
%REGULARIZE Collapse the segments shorter than a given threshold into vertices.
%	OBJ = REGULARIZE(OBJ, EPSILON) removes every segment whose curvilinear
%	length is lower than EPSILON. Each removed segment is replaced by a
%	single vertex located at the barycenter of its knots; the neighbouring
%	segments, the grain loops and the interfaces are updated accordingly.
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

	n_new = length(short_segment);
	new_vertices = zeros(n_new, size(V, 2));
	for i = 1:n_new
		new_vertices(i, :) = mean(V(seg{short_segment(i)}, :), 1);
	end
	new_vertices_ids = (1:n_new) + size(V, 1);

	grains = obj.Grains;
	outerLoops = grains.OuterLoop;
	innerLoops = grains.InnerLoops;
	for i = 1:height(grains)
		[outerLoops{i}, seg] = regularizeLoop(outerLoops{i}, short_segment, new_vertices_ids, seg);
		inner = innerLoops{i};
		for j = 1:length(inner)
			[inner{j}, seg] = regularizeLoop(inner{j}, short_segment, new_vertices_ids, seg);
		end
		innerLoops{i} = inner;
	end

	% Remove the short segments and remap all the IDs
	keep = true(n_seg, 1);
	keep(short_segment) = false;
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

function [loop, seg] = regularizeLoop(loop, short_segment, new_vertices_ids, seg)
%REGULARIZELOOP Replace each occurrence of the short segments in LOOP by a
%	single vertex spliced into the neighbouring segments.

	k = 1;
	while k <= length(short_segment)
		a = find(abs(loop) == short_segment(k), 1, 'first');
		if isempty(a) || length(loop) <= 2
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
		id_new = new_vertices_ids(k);
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
