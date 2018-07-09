function Ife=FEmisfituv(ures,vres,coordinates,connectivity,nip)
	
	
	% Consider doing this as follows
	% M=MassMatrixBlockDiagonal2D(coordinates,connectivity,nip,CtrlVar);
	% [u;v]'*M*[u;v]/2
    % using M as an input
	
   uCost=FE_inner_product(ures,ures,coordinates,connectivity,nip);
   vCost=FE_inner_product(vres,vres,coordinates,connectivity,nip);

    
    Ife=(uCost+vCost)/2;
    
    
end



