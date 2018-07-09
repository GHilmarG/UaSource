function h = hfunforwardC(x,y,MeshSize)
	
	% User defined size function for square
	
	%h = 0.01 + 0.1*sqrt( (x-0.25).^2+(y-0.75).^2 );
	%h =  2.5e3 +1.0*sqrt(x.*x+y.*y);
	

	 % hgrad=0.0; vgrad=0.2; h =  MeshSize+sqrt(x.*x+y.*y)*hgrad+vgrad*abs(y); h(h>MeshSize*10)=MeshSize*10;
	 h =  zeros(numel(x),1)+MeshSize ;
	
	
end      % hfun1()

