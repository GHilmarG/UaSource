function [IntEle]=fIntElementwise(connectivity,coordinates,nip,f,CtrlVar)
    
  
    %
    % calculates int f where f_p=f_q N_q  F_p with F_p=1 for ele p, zero
    % otherwise, with p=1...Nele
    %
    
    
    
    [Nele,nod]=size(connectivity) ; ndim=2;
    fnod=reshape(f(connectivity,1),Nele,nod);   % Nele x nod
    [points,weights]=sample('triangle',nip,ndim);
    
    
    IntEle=zeros(Nele,1); EleArea=zeros(Nele,1);
    
    for Iint=1:nip
        
        
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        [~,detJ]=derivVector(coordinates,connectivity,nip,Iint);
        fint=fnod*fun;
        
        detJw=detJ*weights(Iint);
        EleArea=EleArea+detJw;
        
        IntEle=IntEle+fint.*detJw;
        
    end
    
end





