function n = grainNumPixel(grains)
%%GRAAINNUMPIXEL Number of pixels (measurements) in each grain.
%
%	n = GRAAINNUMPIXEL(grains) is compatible with both MTEX 5.x and 6.0,
%	where the number of pixels is stored in grains.grainSize, and MTEX 6.1
%	or later, where it has been renamed to grains.numPixel.
%
%	See also grainVertices.
	if isprop(grains,'numPixel')	% MTEX >= 6.1
		n = grains.numPixel;
	else							% MTEX <= 6.0
		n = grains.grainSize;
	end
end
