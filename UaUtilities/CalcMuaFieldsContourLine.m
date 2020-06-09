function [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,Field,Value,varargin)

    
    if isempty(Field)  || numel(Field)~=MUA.Nnodes
        xc=[] ; yc=[] ; 
        return
    end

    
   
    GF.node=Field ;
    
    CtrlVar.GLthreshold=Value; 
    CtrlVar.GLsubdivide=1 ; 
    CtrlVar.PlotGLs=false ;
    
    
    [xc,yc]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],varargin{:}) ;
    
    
    
    
end