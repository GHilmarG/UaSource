
function BCs=GetBoundaryConditions(Experiment,CtrlVar,MUA,BCs,time,s,b,h,S,B,ub,vb,ud,vd,GF)

narginchk(15,15)

% does DefineBoundaryConditions.m exist in the run directory?
% if so then use that instead of DefineBCs.m


if exist(fullfile(cd,'DefineBoundaryConditions.m'),'file')

    fprintf(' Using DefineBoundaryConditions.m to define boundary conditions \n')
    BCs=DefineBoundaryConditions(Experiment,CtrlVar,MUA,BCs,time,s,b,h,S,B,ub,vb,ud,vd,GF);
    
    
else
    
    fprintf(' Using DefineBCs.m to define boundary conditions \n')
    [ubFixedNode,ubFixedValue,vbFixedNode,vbFixedValue,...
        ubTiedNodeA,ubTiedNodeB,vbTiedNodeA,vbTiedNodeB,...
        hFixedNode,hFixedValue,hTiedNodeA,hTiedNodeB]=...
        DefineBCs(Experiment,CtrlVar,MUA,time,s,b,h,S,B,ub,vb,ud,vd,GF);
    
    BCs.ubFixedNode=ubFixedNode;
    BCs.ubFixedValue=ubFixedValue;
    BCs.vbFixedNode=vbFixedNode;
    BCs.vbFixedValue=vbFixedValue;
    
    BCs.ubTiedNodeA=ubTiedNodeA;
    BCs.ubTiedNodeB=ubTiedNodeB;
    BCs.vbTiedNodeA=vbTiedNodeA;
    BCs.vbTiedNodeB=vbTiedNodeB;
    
    BCs.hFixedNode=hFixedNode;
    BCs.hFixedValue=hFixedValue;
    BCs.hTiedNodeA=hTiedNodeA;
    BCs.hTiedNodeB=hTiedNodeB;
    
end

switch lower(CtrlVar.FlowApproximation)
    
    case 'sstream'
        
        BCs.udFixedNode=[];
        BCs.udFixedValue=[];
        BCs.vdFixedNode=[];
        BCs.vdFixedValue=[];
        BCs.udTiedNodeA=[];
        BCs.udTiedNodeB=[];
        BCs.vdTiedNodeA=[];
        BCs.vdTiedNodeB=[];
        BCs.udvdFixedNormalNode=[];
        BCs.udvdFixedNormalValue=[];
        
    case 'ssheet'
        
        BCs.ubFixedNode=[];
        BCs.ubFixedValue=[];
        BCs.vbFixedNode=[];
        BCs.vbFixedValue=[];
        BCs.ubTiedNodeA=[];
        BCs.ubTiedNodeB=[];
        BCs.vbTiedNodeA=[];
        BCs.vbTiedNodeB=[];
        BCs.ubvbFixedNormalNode=[];
        BCs.ubvbFixedNormalValue=[];
        
end


% do some basic checks

if numel(BCs.ubFixedNode) ~=  numel(BCs.ubFixedValue)
    error('GetBoundaryConditions:ub','Number of fixed ub nodes not equal to number of ub fixed values! \n')
end


if numel(BCs.vbFixedNode) ~=  numel(BCs.vbFixedValue)
    error('GetBoundaryConditions:vb','Number of fixed vb nodes not equal to number of vb fixed values! \n')
end



end





