function [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,Field,Value,varargin)

    
    %%
    %
    %   [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,Field,Value,varargin)
    %
    %
    % Calculates the contour line of the nodal variable 'Field' for the value 'Value'
    %
    %  
    %
    %
    % Example: Calculate, and then plot, the ice-thickness contour line for the ice-thickness value 100 
    %
    %
    %   [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F.h,100)
    %   plot(xc,yc) ; axis equal
    %%
    
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