function [detEle]=EleDeterminant(coordinates,connectivity,nip)

    Nele=size(connectivity,1) ;
    
    detEle=zeros(Nele,nip);
    
    for Iint=1:nip
        [~,d]=derivVector(coordinates,connectivity,nip,Iint);
        detEle(:,Iint)=d;
    end
    
end



