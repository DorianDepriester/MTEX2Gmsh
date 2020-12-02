function new_paths=EulerPath(F,datatype)
%EULERPATH Compute the Euler path of a graph.
%
%	EULERPATH(F,DTYPE) where F is an N-2 array defining the connectivity
%	between nodes (defined by their indices). The return value is a cell 
%	array defined by nodes indices written in the DTYPE datatype (can be
%	'uint8', 'uint16', 'uint32' or 'uint64').
%
%	See also computeSegments


%% Take advantage of MTEX's EulerCycles function
paths=EulerCycles(F)';

%% Get rid of ambiguous paths (self-touching boundaries)
% Here, we want to find the shortest path, whereas EulerCycles gives the
% longest one.
new_paths=paths;
for k=1:length(paths)
	long_path=paths{k};
	short_path=long_path(1:end-1);
	u_path = unique(short_path);
	dupl_id = find(hist(short_path, u_path) > 1);	% Looking for duplicates
	n_dupl=length(dupl_id);
	if n_dupl
		new_loops = cell(n_dupl,1);
		
		% We need to fix small loops first, so we may track them first
		dupl_dist = zeros(n_dupl,1);	% Distance between duplicates
		for i=1:n_dupl
			dupl_dist(i)=diff(find(short_path==u_path(dupl_id(i))));
		end
		[~, order] = sort(dupl_dist);	% Deal with small loops first
		for i=1:n_dupl
			dupl_val = u_path(dupl_id(order(i)));	% Duplicate node number
			n_vtx= length(short_path);
			shortcut = [find(short_path==dupl_val,1) find(short_path==dupl_val,1, 'last')];
			
			% Move the loop from initial path to new ones
			new_loops{i} = cast(short_path(shortcut(1): shortcut(2))',datatype);
			short_path=short_path([1:shortcut(1) (shortcut(2)+1):n_vtx]);
		end
		
		% Update the current path and append the new ones
		new_paths{k} = cast([short_path long_path(end)]',datatype);
		new_paths=vertcat(new_paths, new_loops);
	else
		new_paths{k}=cast(new_paths{k}',datatype);
	end
end

%% Sort loops such that the (real) outer loop appears first
loop_len=cellfun('length',new_paths);
[~,order]=sort(loop_len,'descend');	% The largest one is considered as the outer loop
new_paths=new_paths(order);	
