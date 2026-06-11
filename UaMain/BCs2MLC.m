function MLC=BCs2MLC(CtrlVar,MUA,BCs)

narginchk(3,3)



MLC=MultiLinearConstraints;

if isempty(BCs)
    return
end

%% input tests

if numel(BCs.ubTiedNodeA) ~= numel(BCs.ubTiedNodeB) ; save TestSave ; error(' number of elements  in BCs.uTiedNodeA and BCs.uTiedNodeB  not the same ') ; end
if numel(BCs.vbTiedNodeA) ~= numel(BCs.vbTiedNodeB) ; save TestSave ; error(' number of elements  in BCs.uTiedNodeA and BCs.uTiedNodeB not the same ') ; end

% get rid of duplicate boundary conditions and just ignore extra BCs
[ubFixedNodet,itemp]=unique(BCs.ubFixedNode) ; BCs.ubFixedValue=BCs.ubFixedValue(itemp);
[vbFixedNodet,itemp]=unique(BCs.vbFixedNode) ; BCs.vbFixedValue=BCs.vbFixedValue(itemp);
%
% if numel(ubFixedNodet) ~= numel(BCs.ubFixedNode)  ; disp(' Duplicate Dirichlet BCs for u') ; end
% if numel(vbFixedNodet) ~= numel(BCs.vbFixedNode)  ; disp(' Duplicate Dirichlet BCs for v') ; end

BCs.ubFixedNode=ubFixedNodet; BCs.vbFixedNode=vbFixedNodet;

%% now set up Lagrange matrices
% velocity
[ubvbL,ubvbRhs]=CreateLuv(MUA,BCs.ubFixedNode,BCs.ubFixedValue,BCs.vbFixedNode,BCs.vbFixedValue,BCs.ubTiedNodeA,BCs.ubTiedNodeB,BCs.vbTiedNodeA,BCs.vbTiedNodeB,BCs.ubvbFixedNormalNode,BCs.ubvbFixedNormalValue);
[udvdL,udvdRhs]=CreateLuv(MUA,BCs.udFixedNode,BCs.udFixedValue,BCs.vdFixedNode,BCs.vdFixedValue,BCs.udTiedNodeA,BCs.udTiedNodeB,BCs.vdTiedNodeA,BCs.vdTiedNodeB,BCs.udvdFixedNormalNode,BCs.udvdFixedNormalValue);

% thickness
% here add pos. thickness constraints

[hL,hRhs]=createLh(MUA.Nnodes,[BCs.hFixedNode;BCs.hPosNode],[BCs.hFixedValue;BCs.hPosValue],BCs.hTiedNodeA,BCs.hTiedNodeB);

[LSFL,LSFRhs]=createLh(MUA.Nnodes,BCs.LSFFixedNode,BCs.LSFFixedValue,BCs.LSFTiedNodeA,BCs.LSFTiedNodeB);


MLC.ubvbL=ubvbL ; MLC.ubvbRhs=ubvbRhs ;
MLC.udvdL=udvdL ; MLC.udvdRhs=udvdRhs ;
MLC.hL=hL; MLC.hRhs=hRhs;
MLC.LSFL=LSFL; MLC.LSFRhs=LSFRhs;

%LastBCs=BCs ; LastMLC=MLC;

%% scale L

[MLC.ubvbL,MLC.ubvbRhs,isLLubvb]=ScaleL(CtrlVar,MLC.ubvbL,MLC.ubvbRhs) ; 
[MLC.udvdL,MLC.udvdRhs,isLLudvd]=ScaleL(CtrlVar,MLC.udvdL,MLC.udvdRhs) ; 
[MLC.hL,MLC.hRhs,isLLh]=ScaleL(CtrlVar,MLC.hL,MLC.hRhs) ; 
[MLC.LSFL,MLC.LSFRhs,isLLLSF]=ScaleL(CtrlVar,MLC.LSFL,MLC.LSFRhs) ; 


end


