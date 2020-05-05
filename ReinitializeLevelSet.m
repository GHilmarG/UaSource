function [LSF,UserVar,RunInfo]=ReinitializeLevelSet(UserVar,RunInfo,CtrlVar,MUA,LSF)
    
    
    %%
    %   [LSF,UserVar,RunInfo]=ReinitializeLevelSet(UserVar,RunInfo,CtrlVar,MUA,F,LSF)
    %
    %
    % Reinitilizes the Level Set;
    %
    % On return LSF is based on the distance from those nodes that were 'On' the zero-line of
    % the level set as given as input.
    %
    % This distance is signed, and the sign is positive over 'In' nodes and negative over
    % 'Out' nodes.  Hence, this is only usefull for REinitilizing the level set, but can not use to
    % initilize it.
    %
    % In, Out, and On nodes as based on a nodal mask of the Level Set given as input.
    %
    % The value of the On nodes is not changed.
    %
    %%
    
    %%
    % 1) mask
    Mask=CalcMeshMask(CtrlVar,MUA,LSF,0);
    
    % 2) Distance from nodes to the nodes of the zero-line elements
    if numel(find(Mask.NodesOn))>0
        Dist=pdist2(MUA.coordinates(Mask.NodesOn,:),MUA.coordinates,'euclidean','Smallest',1) ;
        
  
        % 3) Replace LSF with signed distance over In and Out nodes
        LSF(Mask.NodesIn)=Dist(Mask.NodesIn) ;
        LSF(Mask.NodesOut)=-Dist(Mask.NodesOut) ;
    else
        
        fprintf('ReinitializeLevelSet:No calving-front nodes found within the domain.\n')
        
    end
    
    
    
end