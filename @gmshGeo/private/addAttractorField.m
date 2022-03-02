function id_field = addAttractorField(ffid, id_field, n_segments, def_size, slope)
	fprintf(ffid,'Field[%i] = Distance;\n',id_field);
	fprintf(ffid,'Field[%i].EdgesList ={1:%i};\n',id_field,n_segments);
	id_field=id_field+1;
	fprintf(ffid,'Field[%i] = MathEval;\n',id_field);
	fprintf(ffid,'Field[%i].F = "F%i*%g+%g";\n',id_field,id_field-1,slope,def_size);
end

