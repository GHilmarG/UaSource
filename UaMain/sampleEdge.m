function [s,wt] = sampleEdge(element,nip,ndim,iEdge)

	% returns sample points and weigths for a line element defined for 0 to 1
	
	%
	% does this by taking sample and weights for a line element defined along -1 to 1 and 
	% shifting it into 0 to 1 and dividing weights by 2
	% only usefull for line elements (good for edges)
	ndim=ndim-1;
	[sEdge,wt] = sample(element,nip,ndim);
	 
	switch element
		case 'line'
			sEdge=(sEdge+1)/2;
			wt=wt/2;
		otherwise
			error('not implemented')
	end
	
	s=zeros(nip,2);
	switch iEdge
		case 1
			s(:,1)=1-sEdge; 
		case 2
			s(:,2)=sEdge; 
		case 3
			s(:,1)=sEdge; s(:,2)=1-sEdge;
		otherwise
			error(' what edge? ')
	end
	
	
end

