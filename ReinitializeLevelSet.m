
function [LSF,UserVar,RunInfo]=ReinitializeLevelSet(UserVar,RunInfo,CtrlVar,MUA,LSF,Threshold,xc,yc)
    
    
    %%
    %   [LSF,UserVar,RunInfo]=ReinitializeLevelSet(UserVar,RunInfo,CtrlVar,MUA,F,LSF)
    %
    %
    % Reinitilizes the Level Set;
    %
    %%

    
    %% 1)
    
    if nargin < 6  || isempty(Threshold)
        Threshold=0 ;
    end
    
    if nargin <7 || (isempty(xc) || isempty(yc)) 
        CtrlVar.LineUpGLs=false ;
        [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,LSF,Threshold);
    end
    
    % 2) Distance 
    if numel(xc)>0
  
        Dist=pdist2([xc(:) yc(:)],MUA.coordinates,'euclidean','Smallest',1) ;
        Dist=Dist(:) ;
    
        PM=sign(LSF) ; 
        LSF=PM.*Dist;
  
        
        
    else
        
        fprintf('ReinitializeLevelSet:No calving-front nodes found within the domain.\n')
        
    end
    
    
    
end