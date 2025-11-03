function jg=FindNearestGroundedPoint(CtrlVar,MUA,s,b,S,B,GF,jf)

%  
% Finds the closest grounded node to nodes in the list jf having smaller draft
%
% For any element i in jf and jf, jg(i) is the closest grounded node to jf(i) 
% such that  b(jg(i)) < b(jf(i))
%
% On return jg is a list of nodes equal in size to jf.
%
% The nodes b(jg) are all grounded.
%

if nargin <8
    jf=find(GF.node<0.5) ;  % These nodes are afloat
end

x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ;



xf=x(jf) ; yf=y(jf) ; % (x,y) coordinates of all floating nodes

% find grounding line, and all elements crossing the grounding line
[GLgeo,~,GLnodes]=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);


xg=x(GLnodes) ; yg=y(GLnodes); % These are locations of nodes belonging to edges 
                               % crossing the grounding line. 

% find distance from all floating nodes to all grounding line nodes
numel(xf)
numel(xg)
[daMin,ia,dbMin,ib,D2]=DistanceCloudAtoCloudB(xf,yf,xg,yg);

% if I only wanted to find nearest grounded node then I could now set 
% jg=GLnodes(ia)
% But I must find a grounded node with greater draft so.
% There is a trade-off between how close a grounded not is to a given floting node
% and much larger the draft of that grounded node than that of the floating node.
% This trade-off is determined by the value of `weigth' parameter.
% Potentially the nearest grounded node with larger
% draft is unreasonably far away, so I allow for the possibility 
% to accept a node with less draft. 
% D2 -> D2+ b(jf)+b(GLnodes)

weight=10; % this affect the relative importance of horizontal distance
            % and difference in draft at the floating and the grounded point.

DraftPenalty=100;  % Increasing this value makes accepting grounded node with less draft
                   % than the floating node, less likely.
                   
[X,Y]=ndgrid(b(jf),b(GLnodes));
DraftDist=Y-X;  
I=DraftDist>0;
DraftDist(I)=DraftPenalty*DraftDist(I); % this is a `signed' distance from the
                                        % draft of the floating node to that of the
                                        % grounded node. It is negative if
                                        % the draft of the grounded node is larger. 
                                        % Since I then find the min distance, this makes
                                        % taking such nodes more likely.
D=sqrt(D2);
xybDist=D + weight*DraftDist;
[xybMinDist2,I]=min(xybDist,[],2);
jg=GLnodes(I);

% in those cases where draft of a floating node is less (greater negative)
% than the draft of the corresponding grounding node, replace
% grounding node with the floating node itslef. This implied no
% plume at that floating point because the thickness of the plume will the be zero
J=b(jf)< b(jg);
jg(J)=jf(J);

% for J=1:numel(jf)
%     if b(jf(J))< b(jg(J))
%         jg(J)=jf(J);
%     end
% end


return
figure
%x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ;
for J=1:1:numel(jf)
    if b(jf(J)) > b(jg(J))
        col='b';
    else
        col='r';
    end
    plot([x(jf(J)) x(jg(J))]/CtrlVar.PlotXYscale ,[y(jf(J)) y(jg(J))]/CtrlVar.PlotXYscale,...
        col,'LineWidth',1)
    hold on
end

%plot(x(GLnodes)/CtrlVar.PlotXYscale,y(GLnodes)/CtrlVar.PlotXYscale,'+c')
plot(x(jg)/CtrlVar.PlotXYscale,y(jg)/CtrlVar.PlotXYscale,'or')
plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'k','LineWidth',1);
axis equal
title(sprintf('weight %f ',weight))



end