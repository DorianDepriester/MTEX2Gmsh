classdef LocalSize
	%LOCALSIZE helps defining local sizes and local size gradients in
	%specific grains. The constructed object is intended to be passed to
	%the 'LocalSize' argument of the MESH command from gmshGeo class.
	%
	% See also gmshGeo/mesh	
	properties
		grainID=[];
		sizeAtBoundaries=[];
		slope=[];
	end
	
	methods
		function obj=LocalSize(id,def_size,varargin)
			%LOCALSIZE Object contructor for class LocalSize.
			% LOCALSIZE(A, S) sets the element size to S in grains with IDs
			% A. S can be a single value or an array of the same size as A.
			%
			% LOCALSIZE(A, S, q) sets the size at grain boundaries to S and 
			% prescribes a size gradient with slope q. Again, q can be a
			% single value of an array of the same size as A.
			%
			%
			% Examples:
			%	ls=LOCALSIZE(5, 0.1) sets the element size to 0.1 in grain	
			% labelled 5.
			%
			%	ls=LOCALSIZE([5 10], 0.1) sets the element size to 0.1 in 
			% grain	labelled 5 and 10.
			%
			%	ls=LOCALSIZE([5 10], [0.1 0.2]) sets the element size to
			%	0.1 in grain labelled 5, and 0.2 in grain labelled 10.
			%
			%	ls=LOCALSIZE(5, 0.1, 0.5) sets the element in grain 5 to 
			% follow a size gradient with increasing distance from grain 
			% boundary with slope 0.5, starting from 0.1.
			%
			%	ls=LOCALSIZE([5 10], 0.1, [0.5 0.7]) does the same as
			%	before, plus it sets the slope in grain labelled 10 to 0.7.
			%
			% See also gmshGeo/mesh
			if nargin
				n=length(id);
				if numel(def_size)==1
					def_size=repmat(def_size,n,1);
				elseif length(def_size)~=n
					error('The sizes at boundaries must be either a unique value or an array of the same length as the grain IDs.')
				end
				if nargin==3
					slope=varargin{1};				
				else
					slope=0.0;
				end
				if numel(slope)==1
					slope=repmat(slope,n,1);
				elseif length(slope)~=n
					error('The size of the slopes must be either a unique value or an array of the same length as the grain IDs.')
				end				
				obj(n)=obj;
				for i=1:n
					obj(i).grainID=id(i);
					obj(i).sizeAtBoundaries=def_size(i);		
					obj(i).slope=slope(i);
				end
			end
		end
		
		function disp(obj)
			t=table([obj(:).grainID]', [obj(:).sizeAtBoundaries]', [obj(:).slope]', 'VariableNames', {'GrainID', 'SizeAtBoundaries', 'Slope'});
			if isempty(t)
				disp('Empty set for local element sizes')
			else
				disp(t);
			end
		end
	end	
end

