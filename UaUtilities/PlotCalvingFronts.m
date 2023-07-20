function [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F,varargin)
    
    %
    % Plots calving fronts, just a simple wrapper around PlotGroundingLines
    %
    %
    % To suppress plotting set CtrlVar.PlotGLs=false
    %
    %
    % F is the Ua field variable, but it can also be given as LSF or any other nodal variable 
    %
    % Example:
    %
    % figure ; [xC,yC]=PlotCalvingFronts(CtrlVar,MUA,F,color="r");
    %
    % figure ; [xC,yC]=PlotCalvingFronts(CtrlVar,MUA,LSF,color="r");
    %
    %
    %  figure ; [xC,yC]=PlotCalvingFronts(CtrlVar,MUA,LSF,"k--");
    %
    % Also consider using: 
    %
    %   [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,Field,Value,varargin)
   
    if  isnumeric(F) && numel(F)==MUA.Nnodes 
        LSF=F;
    else
        LSF=F.LSF ;
    end
    
    if isempty(LSF)
        xc=[] ; yc=[] ;
        return
    end


    CtrlVar.LineUpGLs=true ;

    GF.node=LSF ;

    CtrlVar.GLthreshold=0;
    [xc,yc]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],varargin{:}) ;

 



end