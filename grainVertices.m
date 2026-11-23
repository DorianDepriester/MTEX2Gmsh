function V = grainVertices(grains)
%%GRAAINVERTICES Full list of the grain vertices, as an N x 2 matrix.
%
%	V = GRAAINVERTICES(grains) returns the coordinates of all the vertices
%	of the grains, indexed the same way as grains.boundary.F, grains.poly
%	and grains.triplePoints.id. It is compatible with both MTEX 5.x, where
%	the vertices are stored in grains.V, and MTEX 6.0 or later, where they
%	are stored as a vector3d in grains.allV (grains.V is not usable there).
%
%	See also setGrainVertices, grainNumPixel, gmshGeo, cond_smooth.
	if isprop(grains,'allV')		% MTEX >= 6.0
        xyz = grains.allV.xyz;
		V = xyz(:,1:2);
	else							% MTEX <= 5.x
		V = grains.V;
	end
end
