mtexdata twins
ebsd=ebsd('indexed');

%% Create random locations outside the RoI
xmin=min(ebsd.x);
xmax=max(ebsd.x);
ymin=min(ebsd.y);
ymax=max(ebsd.y);

med_width=[10 10];  % Width of medium in x and y directions, around the RoI
n_rand=1000;
Xmin=xmin-med_width(1);
Xmax=xmax+med_width(1);
Ymin=ymin-med_width(2);
Ymax=ymax+med_width(2);
X=rand(n_rand,2).*repmat([Xmax-Xmin Ymax-Ymin],n_rand,1)+repmat([Xmin Ymin],n_rand,1);
in_roi=(X(:,1)>xmin & X(:,1)<xmax) & (X(:,2)>ymin & X(:,2)<ymax);
X=X(~in_roi,:);
n_rand_in=size(X,1);

%% Compute ODF and create random orientations
odf=calcDensity(ebsd.orientations);
rand_ori=odf.discreteSample(n_rand_in);

%% Append these orientations to the intial EBSD map
props=ebsd.prop;
fn = fieldnames(props);
for k=1:numel(fn)
    if strcmp(fn{k},'x')
        props.(fn{k})=X(:,1);
    elseif strcmp(fn{k},'y')
        props.(fn{k})=X(:,2);
    else
        props.(fn{k})=NaN(n_rand_in,1);
    end    
end
fake_ebsd=EBSD(rotation(rand_ori),2*ones(n_rand_in,1),ebsd.CS,props,'unitCell',ebsd.unitCell);
merged_ebsd=[ebsd fake_ebsd];

%% Calc grains
grains_with_medium=calcGrains(merged_ebsd);
plot(grains_with_medium,grains_with_medium.meanOrientation)

% Just for comparison, plot the grains when no medium is added
grains_no_medium=calcGrains(ebsd);
figure
plot(grains_no_medium,grains_no_medium.meanOrientation)


