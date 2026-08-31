function MLC=BCs2MLC(CtrlVar,MUA,BCs)


persistent nInfo


if isempty(nInfo)
    nInfo=0;
end

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

%% Optional RowSubsetSelection

if ~isfield(CtrlVar,"BCsRowSubsetSelection")
    CtrlVar.BCsRowSubsetSelection=false;
end

if CtrlVar.BCsRowSubsetSelection

    if nInfo==0
        fprintf("BCs2MLC: Checking if boundary conditions are linearly independent.\n")
    end


    [ubvbL,row_idx,flag] = RowSubsetSelection(ubvbL);
    if flag==1
        ubvbRhs=ubvbRhs(row_idx);
        if nInfo==0
            fprintf("BCs2MLC: Found linear dependent (ub,vb) boundary conditions. This was corrected by selecting an independent subset. \n")
        end
    end

    [udvdL,row_idx,flag] = RowSubsetSelection(udvdL);
    if flag==1
        udvdRhs=udvdRhs(row_idx);
        if nInfo==0
            fprintf("BCs2MLC: Found linear dependent (ud,vd) boundary conditions. This was corrected by selecting an independent subset. \n")
        end
    end

    [hL,row_idx,flag] = RowSubsetSelection(hL);
    if flag==1
        hRhs=hRhs(row_idx);
        if nInfo==0
            fprintf("BCs2MLC: Found linear dependent h boundary conditions. This was corrected by selecting an independent subset. \n")
        end
    end

    [LSFL,row_idx,flag] = RowSubsetSelection(LSFL);
    if flag==1
        LSFRhs=LSFRhs(row_idx);
        if nInfo==0
            fprintf("BCs2MLC: Found linear dependent level-set boundary conditions. This was corrected by selecting an independent subset. \n")
        end
    end

    nInfo=nInfo+1;
    if nInfo==1
        fprintf("BCs2MLC: No further information about correcting for linear dependent BCs will be provided. \n")
        fprintf("BCs2MLC: and any further corrections will be done silently. \n")
    end
end

%%


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


