function [der]=shape_derEdge(Iint,ndim,nod,points)
	
	if ndim ~=1 ; error(' input not valid ' ) ; end
	
	points=2*points-1;
	
	[der]=shape_der(Iint,ndim,nod,points);

	der=2*der;
	
end
	