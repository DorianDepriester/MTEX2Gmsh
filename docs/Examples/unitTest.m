% This file is a proof of concept for the gmshGeo class. It provides a
% couple of examples, taken from the dataset provided along with MTEX.
%
% This script automatically creates the mesh files from the investigated
% geometries. In addition, it plots the B-splne approximation of Grain 
% Boundaries (GBs) ontop of the computed grains, evidencing the high 
% fidelity of such approximation.

close all

examples = {'small', 'aachen', 'titanium', 'twins'};	% Example datasets

for i=1:length(examples)
	%% MTEX stuff
	% Load EBSD data and compute the grains
	eval(['mtexdata ' examples{i}]);	% Load example data into variable 'ebsd'
	ebsd=ebsd('indexed');				% Remove unindex points
	grains=calcGrains(ebsd);			% Compute grains
	
	%% Process the grains
	% Construct gmshGeo objects differents ways depending on the dataset
	switch examples{i}
		case 'small'
			% Use default meshing paremeters (Wedge elements, element size
			% equal to EBSD resolution)
			G=gmshGeo(grains);
			mesh(G,'small.msh')
		
		case 'aachen'
			% Use simplified GBs and apply size gradient to elements
			G=gmshGeo(grains);
			G=simplify(G);
			mesh(G,'aachen.msh','ElementSize',0.7,'gradient',0.2)
		
		case 'titanium'
			% Smooth the GBs and add a medium surrounding the ROI
			G=gmshGeo(cond_smooth(grains));
			mesh(G,'titanium_medium.msh','elementSize',40,'medium',[2000 2000 200],'mediumElementSize',180);
		
		otherwise
			% Smooth the GBs and mesh using brick element and size gradient
			G=gmshGeo(cond_smooth(grains));
			mesh(G,'twins.msh','ElementSize',0.2,'gradient',0.5,'elementType','Brick');	
	end
	
	%% Plot Grain Boundaries for comparison
	% Visual comparison between the GBs given by MTEX and the approximated
	% ones (B-spline approximations).
	h=figure(i);
	if length(grains.mineralList)==2	% single-phase EBSD map
		% Coulour grains depending their IDs
		plot(grains,grains.id,'noBoundary')
		colormap lines
	else
		% Coulour grains depending on their phase
		plot(grains,'noBoundary')
	end
	hold on
	plot(G)
	legend('Location', 'NorthWest')
	set(h,'Name',examples{i})
end