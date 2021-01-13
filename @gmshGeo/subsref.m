function varargout=subsref(obj,s)
   %	obj(i) only selects the data related to the i-th grain
	switch s(1).type
		case '.'
			if nargout
				varargout{:}=builtin('subsref',obj,s);
			else
				builtin('subsref',obj,s);
			end
		case '()'
			sref=obj;
			grain_tab=sref.Grains;
			k = s.subs;
			if length(k)>1
				error('Only single index can be used here. Consider using an array of indices instead.')
			end						

			%%	Select the grains in the table
			if all(cellfun(@(x) isnumeric(x),k))
				rows=k{:};
				if any(rows<1)
					error('Indices should be stricly positive.')
				end
				if any(rows>height(grain_tab))
					error('The index is larger than the number of grains (%i here).',height(grain_tab))
				end
			else
				rows=strcmp(grain_tab.Phase,k);
			end
			grain_tab=grain_tab(rows,:);
			sref.Grains=grain_tab;

			%%	Keep only the related segments
			if isempty(grain_tab)
				sref.Interfaces=[];
			else
				Out=grain_tab{:,3};
				Out_segIDs=abs(vertcat(Out{:})); %	Concatenate loop-wise
				In=grain_tab{:,4};
				In=vertcat(In{:});				%	Concatenate grain-wise
				In_segIDs=abs(vertcat(In{:}));	%	Concatenate loop-wise 
				segIDs=unique([Out_segIDs; In_segIDs]);

				%%	Update the interfaces
				interfaces=sref.Interfaces;
				for i=1:height(interfaces)
					interface=cast(interfaces.SegmentIDs{i},'like',segIDs);
					interfaces.SegmentIDs{i}=interface(ismember(interface,segIDs));
				end
				t=cellfun('isempty',interfaces.SegmentIDs);	% Keep only non empty sets of segments
				sref.Interfaces=interfaces(~t,:);
			end
		case '{}'
			error('gmshGeo:subsref',...
				'Not a supported subscripted reference')
	end
end