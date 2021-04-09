function id_field = addRestrictField(ffid, id_field, surface_list)
	fprintf(ffid,'Field[%i] = Restrict;\n',id_field);
	fprintf(ffid,'Field[%i].InField = %i;\n',id_field,id_field-1);
	str_surface_list=strjoin(strsplit(num2str(surface_list)),',');
	fprintf(ffid,'Field[%i].SurfacesList = {%s};\n',id_field,str_surface_list);
end

