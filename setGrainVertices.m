function grains = setGrainVertices(grains, V)
%%SETGRAINVERTICES Set the grain vertices from an N x 2 matrix.
%
%	grains = SETGRAINVERTICES(grains,V) updates the vertices of the grains
%	so that their x,y coordinates are given by the N x 2 matrix V, indexed
%	the same way as grains.boundary.F. It is compatible with both MTEX 5.x
%	and MTEX 6.0 or later.
%
%	See also grainVertices.
	if isprop(grains,'allV')		% MTEX >= 6.0
		V3 = grains.allV.xyz;
		V3(:,1:2) = V;
		grains.allV = vector3d(V3(:,1),V3(:,2),V3(:,3));
	else							% MTEX <= 5.x
		grains.V = V;
	end
end
