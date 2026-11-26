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
h = figure;
plot(G);

%% Smoothing the grain boundaries
% One can take advantage of the MTEX |smoothBoundary| function before
% constructing the gmshGeo object:
grains_smooth=grains.smoothBoundary;

%%
% The resulting geometry is obviously smoother than the original one.
G_smooth=gmshGeo(grains_smooth);
plot(G_smooth);


%%
% Note that smoothing is also accompanied with the |simplifyBoundary| 
% function, which applies a Douglas-Peucker algorithm to reduce the number 
% of vertices describing the grain boundaries. 
figure
hold on
G_array = {G, G_smooth};
for i=1:2
    n=cellfun('length',G_array{i}.Segments);
    histogram(n)
end
legend('Raw','Smoothed & simplified')
xlabel('Number of vertices per segment')
ylabel('Frequency')

%% Regularization
% Small segments in the constructed geometry can be automatically 
% collapsed. This can be useful if one wants to avoid very small elements 
% when meshing.
G_regul = G_smooth.regularize(20);

%%
% One can check that the regularized geometry is close to the initial one:
figure(h)
hold on
plot(G_regul)

% Still, small edges have indeed been removed:
figure
hold on
G_array = {G, G_smooth, G_regul};
for i=1:2
    L = G_array{i}.segmentLenghts;
    h = histogram(L);
    h.BinWidth = 5;
end
legend('Smooth geometry', 'Smooth & regularized')
xlabel('Segment lengths')
ylabel('Counts')


%% 
% <html><hr></html>
%
% <index.html Go back to documentation index>