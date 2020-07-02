setMTEXpref('generatingHelpMode',true); % Avoid some artefact (fix issue #5)
%%
mtexdata small
ebsd = ebsd('indexed');
grains = calcGrains(ebsd);
G=gmshGeo(grains);

%% Plot the whole geometry
% The whole geometry can be plotted with the usual plot command:
plot(G)
%%
% The orientation used for plotting is inherited from the MTEX preferences.
% It ensures consistency between the two kinds of plot. For instance, if
% you want to superimpose the grains boundaries defined in G and the
% original grains:
plot(grains,'noBoundary')
hold on
plot(G)
legend('Location', 'NorthWest')

%% Plot specific grains
% Linear indexing helps plotting specific grains. For example, the
% following plots grains 65 and 67:
figure
plot(G([65 67]))
legend('Location', 'Best')
%%
% Alternatively, one can specify the phase to be plotted. E.g:
plot(G('Diopside'))
legend('Location', 'Best')

%% 
% <html><hr></html>
%
% <index.html Go back to documentation index>