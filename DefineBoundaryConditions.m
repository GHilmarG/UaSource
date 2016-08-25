function  BCs=DefineBoundaryConditions(UserVar,CtrlVar,MUA,BCs,time,s,b,h,S,B,ub,vb,ud,vd,GF)
%%
% BCs=DefineBoundaryConditions(UserVar,CtrlVar,MUA,BCs,time,s,b,h,S,B,ub,vb,ud,vd,GF)
%
% BC is a matlab object with the following fields 
%
% BCs = 
% 
%   BoundaryConditions with properties:
% 
%              ubFixedNode: []
%             ubFixedValue: []
%              vbFixedNode: []
%             vbFixedValue: []
%              ubTiedNodeA: []
%              ubTiedNodeB: []
%              vbTiedNodeA: []
%              vbTiedNodeB: []
%      ubvbFixedNormalNode: []
%     ubvbFixedNormalValue: []
%              udFixedNode: []
%             udFixedValue: []
%              vdFixedNode: []
%             vdFixedValue: []
%              udTiedNodeA: []
%              udTiedNodeB: []
%              vdTiedNodeA: []
%              vdTiedNodeB: []
%      udvdFixedNormalNode: []
%     udvdFixedNormalValue: []
%               hFixedNode: []
%              hFixedValue: []
%               hTiedNodeA: []
%               hTiedNodeB: []
%                 hPosNode: []
%                hPosValue: []
%       
%
% see also BoundaryConditions.m
%
%
% Example:
% Setting ub at nodes 1 and 3 to zero.
%
%    BCs.ubFixedNode=[ 1 ; 3 ]
%    BCs.ubFixedValue=BCs.ubFixedNode*0;
% 
% 
%
% Example:
% Setting ub and vb along all grounded nodes at the computational boundary to
% zero.
%
%   x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2); 
%   I=find(GF.node(MUA.Boundary.Nodes)>0.5);  
%   BCs.vbFixedNode=MUA.Boundary.Nodes(I);
%   BCs.ubFixedNode=MUA.Boundary.Nodes(I);
%   BCs.ubFixedValue=BCs.ubFixedNode*0;
%   BCs.vbFixedValue=BCs.vbFixedNode*0;
%   figure ;
%   PlotMuaMesh(CtrlVar,MUA);
%   hold on 
%   plot(x(MUA.Boundary.Nodes)/1000,y(MUA.Boundary.Nodes)/1000,'.b')
%   plot(x(MUA.Boundary.Nodes(I))/1000,y(MUA.Boundary.Nodes(I))/1000,'or')   
% 
% Example: Nodal ties. Setting ub and vb at node 10 equal to ub and vb at node 5 
%
%   udTiedNodeA=3;
%   udTiedNodeB:5;
%   vdTiedNodeA=3;
%   vdTiedNodeB=5; 
%
%
%%



% implementing periodic boundary conditions
% find nodes along boundary 
x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
xd=max(x(:)) ; xu=min(x(:)); yl=max(y(:)) ; yr=min(y(:));
nodesd=find(abs(x-xd)<1e-5); [~,ind]=sort(MUA.coordinates(nodesd,2)); nodesd=nodesd(ind);
nodesu=find(abs(x-xu)<1e-5); [~,ind]=sort(MUA.coordinates(nodesu,2)); nodesu=nodesu(ind);
nodesl=find(abs(y-yl)<1e-5); [~,ind]=sort(MUA.coordinates(nodesl,1)); nodesl=nodesl(ind);
nodesr=find(abs(y-yr)<1e-5); [~,ind]=sort(MUA.coordinates(nodesr,1)); nodesr=nodesr(ind);

BCs.ubTiedNodeA=[nodesu;nodesl]; BCs.ubTiedNodeB=[nodesd;nodesr]; 
BCs.vbTiedNodeA=[nodesu;nodesl]; BCs.vbTiedNodeB=[nodesd;nodesr]; 


end