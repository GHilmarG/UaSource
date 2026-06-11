function [M] = connectivity2adjacency(connectivity)

%  creates the (nodal) adjacency matrix
% (i,j)>=1 if node j belongs to the same element as node i (i,j)=N implies that node i is connected to node j in a total of N
% elements
%
%
%   Note: As of R2015b matlab supports graphs, so an adjacency matrix can be now calculated as:
%
%     g = graph(connectivity, connectivity(:, [2 3 1])); M = adjacency(g);
%
% Three is a difference, as this implementation counts how often nodes are connected, and considers nodes to be connected to
% themselves. These can be made identical as well, but this additional information can be useful.
%
% Turns out that the approach here is about equally fast as the graph theory approach.
%
%%

Nnod=max(connectivity(:));
[Nele,nod]=size(connectivity);


N=Nele*nod*nod;
i=zeros(N,1) ; j=i ; s=ones(N,1); nt=0;

for J=1:nod
    for K=1:nod

        nstak=nt*Nele+(1:Nele)';
        i(nstak)=connectivity(1:Nele,J);
        j(nstak)=connectivity(1:Nele,K);
        nt=nt+1;
    end
end

M=sparseUA(i,j,s,Nnod,Nnod);



%% same as graph theory
% diag_M = diag(M);
% MD = sparse(1:length(diag_M), 1:length(diag_M), diag_M);
% M=sparseUA(i,j,s,Nnod,Nnod);
% M = M - MD;
% [i,j,s]=find(M);
% s=s*0+1;
% M=sparseUA(i,j,s,Nnod,Nnod);
%%

% g = graph(connectivity, connectivity(:, [2 3 1]));    M2 = adjacency(g);




end


