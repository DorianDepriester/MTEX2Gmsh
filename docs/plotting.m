%%
%   mtexdata small
%   ebsd = ebsd('indexed');
%   grains = calcGrains(ebsd);
%   G=gmshGeo(grains);

%% Plot the whole geometry
% The whole geometry can be plotted with the usual plot command:
plot(G)

%% Plot specific grains
% Linear indexing helps plotting specific grains. For example, the
% following plots grains 65 and 67:
plot(G([65 67]))
legend('Location', 'Best')