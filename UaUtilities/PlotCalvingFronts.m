function [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F,varargin)
    
    %
    % Plots calving fronts, just a simple wrapper around PlotGroundingLines
    %
    %
    % To suppress plotting set CtrlVar.PlotGLs=false ;
    %
    %
    %
    
    if isempty(F.LSF)
        xc=[] ; yc=[] ; 
        return
    end
    
    GF.node=F.LSF ;
    
    
    [xc,yc]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],varargin{:}) ;
    
    
    
    
end