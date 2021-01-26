function [ grains ] = cond_smooth(grains,varargin)
%COND_SMOOTH Conditional smoothing of the grain boundaries.
%
%	COND_SMOOTH(grains,SP) smoothes the grain boundaries, but it keeps the
%	singular points unchanged. It is based on the smooth method from MTEX;
%	therefore it accepts the same input arguments.
%
%	See also grain2d/smooth singularPoints.
	sp = singularPoints(grains);
	sp = vertcat(sp{:});
	V_old=grains.boundary.V;
    grains=smooth(grains,varargin{:});
    grains.boundary.V(sp,:)=V_old(sp,:);
end

