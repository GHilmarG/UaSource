function z=Hi(x,y,x0,y0,sx,sy,sz)
	
	% (x0,y0) : coordinate shift
	% sx , sy, sz  : coordinate scaling
	% x=(x-x0)/sx; y=(y-y0)/sy;
	
	
	
	x=(x-x0)/sx; y=(y-y0)/sy;
	z = exp(-x.^2-y.^2/2).*cos(4*x) + exp(-3*((x+0.5).^2+y.^2/2));
	z(z>0.001)=0.001;
	z(z<0)=0;
	z=z/0.001;
	z=z/sz;
	
end