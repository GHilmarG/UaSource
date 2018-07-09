function M = AdjacentNodes(connectivity,EleList)

%%
%  M = AdjacentNodes(connectivity,EleList)
%  creates the (nodal) neighbour matrix.
% (i,j)>=1 if node j is a neighbouring node to node i.
% (i,j)=N implies that node j is a neighbouring node to node i in a total of N
% elements.
%
% EleList is an optional argument that can be used to resrict the list to
% selected elements.
% If EleList is not given, then all elements are used.
%
% Example: To find the neighbouring nodes to node nr. 10: 
%
%    M = AdjacentNodes(connectivity)
%    find(M(10,:))
%
%%

Nnod=max(connectivity(:));
[Nele,nod]=size(connectivity);
switch nod
    case 3
        S =[3  1  2  ; ...
            2  3  1  ];
    case 6
        S =[6  1  2  3  4  5 ; ...
            2  3  4  5  6  1];
    case 10
        S =[9  1  2  3  4  5  6  7  8; ...
            2  3  4  5  6  7  8  9  1];
end

if nargin < 2
    EleList=1:Nele;
end

Nele=numel(EleList);

if nod== 10 ; nod=9 ; end  % for the 10 node element the 10th node is never connected

N=Nele*nod*2;
i=zeros(N,1) ; j=i ;s=ones(N,1); nt=0;


for J=1:nod
    for K=1:2

        nstak=nt*Nele+[1:Nele]';
        i(nstak)=connectivity(EleList,J);
        j(nstak)=connectivity(EleList,S(K,J));
        nt=nt+1;
    end
end

M=sparseUA(i,j,s,Nnod,Nnod);

end


