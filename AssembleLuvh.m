
function [L,Lrhs,lambda]=AssembleLuvh(Luv,Lh,Luvrhs,Lhrhs,lambdauv,lambdah,Nnodes)
	
	% Luv  : #uv constrains x 2Nnodes
	% Lh  : #h constrains x Nnodes
	% L=[Luv 0]
	%   [0  Lh]
	
	[nu,mu]=size(Luv) ;
    [nh,mh]=size(Lh) ;
    
	
	
	if isempty(Lh) && ~isempty(Luv)
		L=[Luv sparse(nu,Nnodes)] ; 
		Lrhs=Luvrhs ;
		lambda=lambdauv;
	elseif ~isempty(Lh) && isempty(Luv)
		L=[sparse(nh,2*Nnodes) Lh] ; 
		Lrhs=Lhrhs ;
		lambda=lambdah;
	elseif ~isempty(Lh) && ~isempty(Luv)
       
		L=[ Luv sparse(nu,Nnodes) ; sparse(nh,2*Nnodes) Lh];
		Lrhs=[Luvrhs;Lhrhs];
		lambda=[lambdauv;lambdah];
	else
		L=[] ; Lrhs=[] ; lambda=[];
	end
		
	
end
	
	