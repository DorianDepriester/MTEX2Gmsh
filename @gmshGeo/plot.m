		function varargout=plot(obj,varargin)
		%PLOT Plot the segments found in each grains.
		%	PLOT(obj) plots all the segments.
		%
		%	PLOT(obj(I)) plots the segments of grains whose indices are
		%	given by the array I.
		%
		%	PLOT(obj(P)) plots the segments of grains of phase P only,
		%	where P is a string.
		%
			if isempty(obj.Grains) || isequal(obj.Grains,struct)
				warning('Empty set. Nothing to plot.')
				return
			end
			T=obj.Interfaces;
			intnames=T.Properties.RowNames;
			nInt=length(intnames);
			D=cell(2*nInt,1);
			for i=1:nInt
				lineSet=T{i,1};
				lineSet=vertcat(lineSet{:});
				X=[];Y=[];
				for j=1:length(lineSet)
					lineID=lineSet(j);
					Vids=obj.Segments{lineID,1};
					x=obj.V(Vids,1);
					y=obj.V(Vids,2);
					if Vids(1)==Vids(end)
						XYbs=BSpline([x(1:end-1),y(1:end-1)],'order',4,'periodic',true);
					else
						XYbs=BSpline([x,y],'order',4);
					end
					X=[X; NaN; XYbs(:,1)]; %#ok<AGROW>
					Y=[Y; NaN; XYbs(:,2)]; %#ok<AGROW>
				end
				D{2*i-1}=X;
				D{2*i}=Y;
			end
		    p=plot(D{:},varargin{:});
			set(p,'LineWidth',2);
			for i=1:nInt
				set(p(i),'DisplayName',intnames{i})
			end
			xlabel('x');
			ylabel('y');
			setPlotOrientation
			if nargout
				varargout{:}=p;
			end
		end
		
function setPlotOrientation
	if strcmpi(getMTEXpref('xAxisDirection'),'east') && strcmpi(getMTEXpref('zAxisDirection'),'intoPlane')
		set(gca,'Ydir','reverse')
	elseif strcmpi(getMTEXpref('xAxisDirection'),'north')
		if strcmpi(getMTEXpref('zAxisDirection'),'outOfPlane')
			set(gca,'Ydir','reverse')
		end
		view([90 -90])
	elseif strcmpi(getMTEXpref('xAxisDirection'),'west')
		set(gca,'Xdir','reverse')
		if strcmpi(getMTEXpref('zAxisDirection'),'outOfPlane')
			set(gca,'Ydir','reverse')
		end
	elseif strcmpi(getMTEXpref('xAxisDirection'),'south')
		view([90 -90])
		set(gca,'Xdir','reverse')
		if strcmpi(getMTEXpref('zAxisDirection'),'intoPlane')
			set(gca,'Ydir','reverse')
		end
	end
	axis equal
end		