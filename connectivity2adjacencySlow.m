function [M] = connectivity2adjacencySlow(connectivity)
    
    Nnod=max(connectivity(:));
    [Nele,nod]=size(connectivity);
    
    iv=zeros(nod,1); jv=iv;
    s=ones(nod,1); M=sparse(Nnod,Nnod);
    
    for I=1:Nele
        for J=1:nod
            iv(1:nod,1)=connectivity(I,J);
            jv(1:nod,1)=connectivity(I,:);
            M=M+sparse(iv,jv,s,Nnod,Nnod);
        end
    end
    
    
end

