
function [L,Lrhs,l]=AssembleLuvh(Luv,Lh,cuv,ch,luv,lh,Nnodes)
	
error('AssembleLuvh:obsolete','AssembleLuvh no longer used. Consider using AssembleLuvhSSTREAM.m instead')

	% Luv  : #uv constrains x 2Nnodes
	% Lh  : #h constrains x Nnodes
	% L=[Luv 0]
	%   [0  Lh]
	
	[nu,~]=size(Luv) ;
    [nh,~]=size(Lh) ;
    
	
	
	if isempty(Lh) && ~isempty(Luv)
		L=[Luv sparse(nu,Nnodes)] ; 
		Lrhs=cuv ;
		l=luv;
	elseif ~isempty(Lh) && isempty(Luv)
		L=[sparse(nh,2*Nnodes) Lh] ; 
		Lrhs=ch ;
		l=lh;
	elseif ~isempty(Lh) && ~isempty(Luv)
       
		L=[ Luv sparse(nu,Nnodes) ; sparse(nh,2*Nnodes) Lh];
		Lrhs=[cuv;ch];
		l=[luv;lh];
	else
		L=[] ; Lrhs=[] ; l=[];
	end
		
	
end
	
	