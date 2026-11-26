function BS = bspline_auto(V)
    if all(V(1, :) == V(end, :))
        if length(V) == 4
            order = 3;
        else
            order = 4;
        end
        BS = BSpline(V(1:end-1,:), 'order', order, 'periodic', true);
    else
        BS = BSpline(V, 'order', 4);
    end
end