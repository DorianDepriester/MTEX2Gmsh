function BS = BSpline(knots,varargin)
%BSPLINE computes the B-spline approximation from a set of coordinates.
% BSPLINE(KNOTS) returns the B-spline interpolation between the reference
% points (knots) whose coordinates are given by the array KNOTS.
% The coordinates of the knots are given vertically, i.e. KNOTS(i,j) gives 
% the j-th coordinate of the i-th knot. The knots can be of any dimension.
%
% BSPLINE(KNOTS,'order',n) Uses n -th order approximation (default: n=2)
%
% BSPLINE(KNOTS,'nint',m) gives m points per interval (default: m=10)
%
% If KNOTS is of size [p,q], the result will be of size [(p-1)*(m-1)+1 ,q],
% except if periodicity is requested (see below).
%
% BSPLINE(KNOTS,'periodic',true) use periodic conditions at end knots.
%

	ip = inputParser;
	addOptional(ip,'order',2)
	addOptional(ip,'nint',10)
	addOptional(ip,'periodic',false)
	parse(ip,varargin{:});
	
	if ip.Results.periodic
		np_rep=ip.Results.order+1;
		knots=[knots(end-np_rep+1:end,:); knots; knots(1:np_rep,:)];
	end	
	
	p=size(knots,1);
	q=size(knots,2);
	
 	if p<=2
 		BS=knots;
 		return
	end
	
	n=ip.Results.nint;
	n=(n-1)*(p-1)+1;	% Overall number of queried points
	y = linspace(0,1,n);

	order=min(ip.Results.order+1,p);
	Xl = zeros(order,order,q);
	t = [zeros(1,order-1),linspace(0,1,p-order+2),ones(1,order-1)];

	BS=zeros(n,q);
	m_min=1;
	m_max=n;
	for m = 1:n-1
		t0 = y(m);
		k = find(t0 >= t,1,'last');
		if (k > p)
			BS=BS(1:m-1,:);
			return;
		end
		Xl(:,1,:) = knots(k-order+1:k,:);
		if k<=order+1
			m_min=max(m,m_min);
		end
		if k>=p-order+2
			m_max=min(m,m_max);
		end

		for i = 2:order
			for j = i:order
				num = t0-t(k-order+j);
				if num == 0
					wt = 0;
				else
					s = t(k+j-i+1)-t(k-order+j);
					wt = num/s;
				end
				Xl(j,i,:) = (1-wt)*Xl(j-1,i-1,:) + wt*Xl(j,i-1,:);
			end
		end
		BS(m,:) = Xl(order,order,:);
	end
	BS(end,:)=knots(end,:);

	if ip.Results.periodic
		BS=BS(m_min:m_max-1,:);
		BS(end,:)=BS(1,:);
	end
end


