function s=printMin(a)
	if ischar(a)
		s=a;
	elseif iscell(a)
		if length(a)==1
			s=a{1};
		else
			s=strcat(sprintf('Min(%s,',a{1}), printMin(a(2:end)), ')');
		end
	end
end