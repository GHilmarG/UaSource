
function [L,Lrhs,l]=AssembleLuvh(Luv,Lh,Luvrhs,Lhrhs,luv,lh,Nnodes)
	
	% Luv  : #uv constrains x 2Nnodes
	% Lh  : #h constrains x Nnodes
	% L=[Luv 0]
	%   [0  Lh]
	
	[nu,~]=size(Luv) ;
    [nh,~]=size(Lh) ;
    
	
	
	if isempty(Lh) && ~isempty(Luv)
		L=[Luv sparse(nu,Nnodes)] ; 
		Lrhs=Luvrhs ;
		l=luv;
	elseif ~isempty(Lh) && isempty(Luv)
		L=[sparse(nh,2*Nnodes) Lh] ; 
		Lrhs=Lhrhs ;
		l=lh;
	elseif ~isempty(Lh) && ~isempty(Luv)
       
		L=[ Luv sparse(nu,Nnodes) ; sparse(nh,2*Nnodes) Lh];
		Lrhs=[Luvrhs;Lhrhs];
		l=[luv;lh];
	else
		L=[] ; Lrhs=[] ; l=[];
	end
		
	
end
	
	