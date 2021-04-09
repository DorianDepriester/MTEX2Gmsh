function writeSequence(ffid,title,idx,Seq)
    if ~isempty(Seq)
		if isempty(idx)
		    fprintf(ffid,'%s{%i',title,Seq(1));
		else
			if isnumeric(idx)
				fprintf(ffid,'%s(%i)={%i',title,idx,Seq(1));
			else
				fprintf(ffid,'%s(%s)={%i',title,idx,Seq(1));
			end
		end
		n=length(Seq);
		if n==1
			fprintf(ffid,'};\n');
		else
			pending=[];
			n_pend=0;
			for k=2:n
				if k~=n
					if (Seq(k)==Seq(k-1)+1 || Seq(k)==Seq(k-1)-1) 
						pending=Seq(k);
						n_pend=n_pend+1;
					else
						if isempty(pending)
							fprintf(ffid,',%i',Seq(k));
						elseif n_pend>1
							fprintf(ffid,':%i,%i',pending,Seq(k));
						else
							fprintf(ffid,',%i,%i',pending,Seq(k));
						end
						pending=[];
						n_pend=0;
					end
				else
					if (Seq(k)==Seq(k-1)+1 || Seq(k)==Seq(k-1)-1)
						fprintf(ffid,':%i};\n',Seq(k));
					else
						if isempty(pending)						
							fprintf(ffid,',%i};\n',Seq(k));
						else
							if n_pend>1
								fprintf(ffid,':%i,%i};\n',pending,Seq(k));
							else
								fprintf(ffid,',%i,%i};\n',pending,Seq(k));
							end
						end
					end
				end
			end
		end
    end
end

