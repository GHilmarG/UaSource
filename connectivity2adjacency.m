function [M] = connectivity2adjacency(connectivity)
    
%  creates the (nodal) adjacency matrix
% (i,j)>=1 if node j belongs to the same element as node i 
% (i,j)=N implies that node i is connected to node j in a total of N elements
%
% 
%

    Nnod=max(connectivity(:));
    [Nele,nod]=size(connectivity);
    
    
    N=Nele*nod*nod;
    i=zeros(N,1) ; j=i ;s=ones(N,1); nt=0;
    
    for J=1:nod
        for K=1:nod
            
            nstak=nt*Nele+[1:Nele]';
            i(nstak)=connectivity(1:Nele,J);
            j(nstak)=connectivity(1:Nele,K);
            nt=nt+1;
        end
    end

    M=sparseUA(i,j,s,Nnod,Nnod);
    
end


