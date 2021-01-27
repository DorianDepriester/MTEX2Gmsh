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
	valide_types=vertcat(fieldnames(sp),'all');
	if ~ischar(sp_type) || ~ismember(sp_type,valide_types)
		error('Wrong type for singular points. It can be ''%s'' or ''all''.',strjoin(valide_types,''', '''));
	elseif strcmp(sp_type,'all')
		sp_d=struct2cell(sp);
		sp_d=[sp_d{:}];
	else
		sp_d=sp.(sp_type);
	end
	c=p.Results.Color;
	s=p.Results.Size;
	m=p.Results.Marker;
	h=plot(V(sp_d,1),V(sp_d,2),'linestyle','none','Marker',m,'MarkerFaceColor',c,'MarkerEdgeColor',c,'MarkerSize',s);
	if nargout
		varargout{:}=h;
	end
end