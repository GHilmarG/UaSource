function [F,RunInfo]=ModifyThicknessBasedOnLevelSet(RunInfo,CtrlVar,MUA,F)
    
    % Set ice thickness where the level set is smaller than 1 to minimum thickness
    
    if ~isempty(F.LSF)
        
        Mask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);
        F.h(Mask.NodesOut)=max(CtrlVar.LevelSetMinIceThickness,CtrlVar.ThickMin) ;
        [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);
        
    end
    
end