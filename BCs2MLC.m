function MLC=BCs2MLC(MUA,BCs)

persistent LastMLC LastBCs

if isequal(BCs,LastBCs)
    MLC=LastMLC;
    return
end

MLC=MultiLinearConstraints;

%% input tests

if numel(BCs.ubTiedNodeA) ~= numel(BCs.ubTiedNodeB) ; save TestSave ; error(' number of elements  in BCs.uTiedNodeA and BCs.uTiedNodeB  not the same ') ; end
if numel(BCs.vbTiedNodeA) ~= numel(BCs.vbTiedNodeB) ; save TestSave ; error(' number of elements  in BCs.uTiedNodeA and BCs.uTiedNodeB not the same ') ; end

% get rid of duplicate boundary conditions and just ignore extra BCs
[ubFixedNodet,itemp]=unique(BCs.ubFixedNode) ; BCs.ubFixedValue=BCs.ubFixedValue(itemp);
[vbFixedNodet,itemp]=unique(BCs.vbFixedNode) ; BCs.vbFixedValue=BCs.vbFixedValue(itemp);

if numel(ubFixedNodet) ~= numel(BCs.ubFixedNode)  ; disp(' Duplicated Dirichlet BCs for u') ; end
if numel(vbFixedNodet) ~= numel(BCs.vbFixedNode)  ; disp(' Duplicated Dirichlet BCs for v') ; end

BCs.ubFixedNode=ubFixedNodet; BCs.vbFixedNode=vbFixedNodet;

%% now set up Lagrange matrices
% velocity
[ubvbL,ubvbRhs]=CreateLuv(MUA,BCs.ubFixedNode,BCs.ubFixedValue,BCs.vbFixedNode,BCs.vbFixedValue,BCs.ubTiedNodeA,BCs.ubTiedNodeB,BCs.vbTiedNodeA,BCs.vbTiedNodeB,BCs.ubvbFixedNormalNode,BCs.ubvbFixedNormalValue);
[udvdL,udvdRhs]=CreateLuv(MUA,BCs.udFixedNode,BCs.udFixedValue,BCs.vdFixedNode,BCs.vdFixedValue,BCs.udTiedNodeA,BCs.udTiedNodeB,BCs.vdTiedNodeA,BCs.vdTiedNodeB,BCs.udvdFixedNormalNode,BCs.udvdFixedNormalValue);

% thickness
% here add pos. thickness constraints

[hL,hRhs]=createLh(MUA.Nnodes,[BCs.hFixedNode;BCs.hPosNode],[BCs.hFixedValue;BCs.hPosValue],BCs.hTiedNodeA,BCs.hTiedNodeB);




%% all done, only tests and setting flags left to do

% % set ties flags (not using this anymore when solving the resulting systems, so presumably not needed)
% MLC.ubvbTies=1 ; MLC.udvdTies=1 ; MLC.hTies=1 ;
% if isempty(BCs.ubTiedNodeA) && isempty(BCs.vbTiedNodeB) && isempty(BCs.ubvbFixedNormalNode) ; MLC.ubvbTies=0  ;end
% if isempty(BCs.udTiedNodeA) && isempty(BCs.vdTiedNodeB) && isempty(BCs.udvdFixedNormalNode) ; MLC.udvdTies=0  ;end
% if isempty(BCs.hTiedNodeA) ; MLC.hTies=0 ; end
% 
% isempty(BCs.ubvbFixedNormalNode) 
% 
% 
% if ~MLC.ubvbTies ;
%     [m,n]=size(ubvbL); if ~isequal(ubvbL*ubvbL',sparse(1:m,1:m,1)) ; save TestSave ; error(' ubvbL transpose(ubvbL) expected to be equal to the identity matrix, but is not!') ; end
% end
% 
% if ~MLC.ubvbTies ;
%     [m,n]=size(udvdL); if ~isequal(udvdL*udvdL',sparse(1:m,1:m,1)) ; save TestSave ; error(' udvdL transpose(udvdL) expected to be equal to the identity matrix, but is not!') ; end
% end
% 
% if ~MLC.hTies
%     [m,n]=size(hL); if ~isequal(hL*hL',sparse(1:m,1:m,1)) ; save TestSave  ; error(' hL transpose(hL) expected to be equal to the identity matrix, but is not!') ; end
% end

MLC.ubvbL=ubvbL ; MLC.ubvbRhs=ubvbRhs ;
MLC.udvdL=udvdL ; MLC.udvdRhs=udvdRhs ;
MLC.hL=hL; MLC.hRhs=hRhs; 

LastBCs=BCs ; LastMLC=MLC;


end


