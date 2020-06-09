function [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F,varargin)
    
    %
    % Plots calving fronts, just a simple wrapper around PlotGroundingLines
    %
    %
    % To suppress plotting set 
    %
    %
    %
   
    if  isnumeric(F) && numel(F)==MUA.Nnodes 
        LSF=F;
    else
        LSF=F.LSF ;
    end
    
    if isempty(LSF)
        xc=[] ; yc=[] ; 
        return
    end
    
    
    
   
    GF.node=LSF ;
    
    CtrlVar.GLthreshold=0; 
    [xc,yc]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],varargin{:}) ;
    
    
    
    
end