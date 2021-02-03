function varargout=plotSingularPoints(obj, varargin)
%PLOTSINGULARPOINTS plots the singular points in the geometry.
%	PLOTSINGULARPOINTS(G) plots the singular points as blue dots.
%
%	PLOTSINGULARPOINTS(G,'type',t) plots the singular points of type t. The 
%	type can be:
%		- 'triplePoints'
%		- 'quadruplePoints'
%		- 'corners'
%		- 'doublePointsOnBorder'
%		- 'symmetric'
%	If 'all' is used, all the singular points are plotted (default
%	behaviour).
%
%	PLOTSINGULARPOINTS(G,...,'Marker',M) uses markers of style named M. For 
%	instance, PLOTSINGULARPOINTS(G,...,'Marker','d') uses diamonds.
%
%	PLOTSINGULARPOINTS(G,...,'Color',C) colors the markers to color named
%	C.
%		
%	PLOTSINGULARPOINTS(G,...,'Size',S) sets the size for the markers
%	(default: 6).
%
%	h=PLOTSINGULARPOINTS(...) return handle to the plot.
%
%	see also plot.
	V=obj.V;
	sp=obj.SingularPoints;
	p = inputParser;
	addOptional(p,'type','all');
	addOptional(p,'Marker','o');
	addOptional(p,'Color','b');
	addOptional(p,'Size',6);
	parse(p,varargin{:});
	
	sp_type=p.Results.type;
	if ~ischar(sp_type)
		error('Type for singular points must be a string.');
	else
		if strcmpi(sp_type,'all')
			sp_d=struct2cell(sp);
			sp_d=vertcat(sp_d{:});
			dname='Singular points';
		else
			if strcmp(sp_type,'triplePoints')
				dname='Triple points';
			elseif strcmp(sp_type,'quadruplePoints')
				dname='Quadruple points';
			elseif strcmp(sp_type,'corners')
				dname='Corners of RoI';
			elseif strcmp(sp_type,'doublePointsOnBorder')
				dname='Double points on border of RoI';
			elseif strcmp(sp_type,'symmetric')
				dname='Symmetric points';					
			else
				valid_types=vertcat(fieldnames(sp),'all');
				error('Wrong type for singular points. It can be ''%s'' or ''all''.',strjoin(valid_types,''', '''));
			end
			sp_d=sp.(sp_type);			
		end
	end
	c=p.Results.Color;
	s=p.Results.Size;
	m=p.Results.Marker;
	h=plot(V(sp_d,1),V(sp_d,2),'linestyle','none','Marker',m,'MarkerFaceColor',c,'MarkerEdgeColor',c,'MarkerSize',s);
	if nargout
		varargout{:}=h;
	end
	set(h,'DisplayName',dname);
end