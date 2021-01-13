function s=size(obj)
	%SIZE Number of grains and size of the ROI.
	s.numberOfGrains=height(obj.Grains);
	vtx=obj.V;
	xmin=min(vtx(:,1));
	xmax=max(vtx(:,1));
	ymin=min(vtx(:,2));
	ymax=max(vtx(:,2));
	s.ROI=[xmax-xmin ymax-ymin];
end	