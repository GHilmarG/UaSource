

function [M] = connectivity2adjacency2(connectivity)
    
    Nnod=max(connectivity(:));
    [Nele,nod]=size(connectivity);
    
    N=Nele*nod*nod;
    i=zeros(N,1) ; j=i ;s=ones(N,1); nstak=0;
    for I=1:Nele
        for J=1:nod
            for K=1:nod
                nstak=nstak+1;
                i(nstak)=connectivity(I,J);
                j(nstak)=connectivity(I,K);
            end
        end
    end
    
    M=sparse(i,j,s,Nnod,Nnod);
    
end


