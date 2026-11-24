setMTEXpref('generatingHelpMode',true); % Avoid some artefact (fix issue #5)
%% Load example data
% Load EBSD data consisting in an hexagonal grid, then compute the grains.
mtexdata titanium
grains = calcGrains(ebsd);
plot(grains,grains.meanOrientation)

%% Compute the serrated geometry
% As you can see, the grain boundaries are heavily serrated because of the 
% hexagonal grid used for EBSD indexing. This drawback obviously also apply
% when constructing the gmshGeo object (despite the smoothing property of 
% the BSpline approximation). Indeed:
G=gmshGeo(grains);	
plot(G);

%% Smoothing the grain boundaries
% One can take advantage of the MTEX |smoothBoundary| function before
% constructing the gmshGeo object:
grains_smooth=grains.smoothBoundary;
plot(grains_smooth,grains_smooth.meanOrientation)

%%
% The resulting geometry is obviously smoother than the original one.
G_smooth=gmshGeo(grains_smooth);
plot(G_smooth);


%%
% Note that smoothing is also accompanied with the |simplifyBoundary| 
% function, which applies a Douglas-Peucker to reduce the number of 
% vertices used for describing the geometry. 
figure
n_o=cellfun('length',G.Segments);
n_s=cellfun('length',G_smooth.Segments);
hist([n_o n_s])
legend({'Raw','Smoothed'})
xlabel('Number of vertices per segment')
ylabel('Frequency')

%% 
% <html><hr></html>
%
% <index.html Go back to documentation index>