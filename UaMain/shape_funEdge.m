function [fun]=shape_funEdge(Iint,ndim,nod,points)
	
	if ndim ~=1 ; error(' input not valid ' ) ; end
	
	points=2*points-1;
	
	[fun]=shape_fun(Iint,ndim,nod,points);

	
end
	