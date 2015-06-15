
%function GLgeo=LineUpEleContours(GLgeo);

load GLgeoFile GLgeo connectivity coordinates


% GLgeo(:,1)                        : list of elements with nodes on both side of the grounding line
% GLgeo(:,2)                        : edge number
% GLgeo(:,[3 4])                    : x y coordinates of grounding line position of element GLgeo(:,1) and edge GL(:,2)
% GLgeo(:,[5 6])                    : x y coordinates of grounding line position of element GLgeo(:,1) and edge GL(:,2)+1  (cyclic)
%
% The  grounding line can be plotted using
%  plot(GLgeo(:,[3 4])',GLgeo(:,[5 6])')
%

% for each element create a pointer towards neighboring element

[FEmeshTriRep]=CreateFEmeshTriRep(connectivity,coordinates);
% FEedges=edges(FEmeshTriRep);
% FEneighbors=neighbors(FEmeshTriRep)
% triplot(FEmeshTriRep,'color','k') ;

edge=[ 1 2 ; 2 3 ; 3 1];

N=1;

Ele=GLgeo(N,1);
EleNodes=edge(GLgeo(1,2),:);


EdgeNodes=FEmeshTriRep.Triangulation(GLgeo(N,1),edge(GLgeo(N,2),:)) ;% nodes of one of the two edges of Ele that the contour goes through

temp=edgeAttachments(FEmeshTriRep,EdgeNodes);
temp{1};
ttt=neighbors(FEmeshTriRep,Ele);
NextEle=intersect(temp{1},ttt)   % this is the element that links

% create a key
key=zeros(max(GLgeo(:,1)),1);
key(GLgeo(:,1))=1:length(GLgeo);


N;
GLgeo(N,[3 4])
GLgeo(N,[5 6])

N=key(NextEle);
GLgeo(N,[3 4])
GLgeo(N,[5 6])
