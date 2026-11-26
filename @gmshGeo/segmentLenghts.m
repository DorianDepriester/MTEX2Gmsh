function lenghts = segmentLenghts(obj)
%SEGMENTLENGTHS computes the curvilinear length of each segment
% SEGMENTLENGTHS first approximate the segments with BSplines, then compute
% the curvilinear length of each curve.
%
% See also regularize

    seg = obj.Segments;
    V = obj.V;
    n_seg = length(seg);
    lenghts = zeros(n_seg, 1);
    for i = 1:n_seg
        Vi = V(seg{i}, :);
        bs = bspline_auto(Vi);
        lenghts(i) = sum(vecnorm(diff(bs, 1, 1), 2, 2));
    end
end
