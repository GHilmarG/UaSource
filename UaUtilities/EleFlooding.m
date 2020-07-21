

function IEle=EleFlooding(CtrlVar,MUA,NodeSeed,EleSubset,CreateFigure)
    
    %%
    %
    %   IEle=EleFlooding(CtrlVar,MUA,NodeSeed,EleSubset,CreateFigure)
    %
    % Finds all elements within an element subset (EleSubset) connected to a given
    % node (NodeSeed)
    %
    % Example: Find all elements downstream of the grounding line that are
    % connected to the node furthest from flotation: 
    %
    %   GF=IceSheetIceShelves(CtrlVar,MUA,GF); xCutOff=400e3; 
    %   EleSubset=GF.ElementsDownstreamOfGroundingLines & (MUA.xEle>xCutOff) ;
    %   hf=rhow*(S-B)./rho ; [MinFlotation,NodeSeed]=min((s-b)-hf);
    %   ElementsToBeDeactivated=EleFlooding(CtrlVar,MUA,NodeSeed,EleSubset);
    %
    % Node: Currently only works for 3-node elements, although this can easily be
    % addressed using the same general idea.  
    %
    %%
    
    if ~islogical(EleSubset)
        error('IEle:InvalidInput','EleSubset must be a logical index')
    end
    
    if islogical(NodeSeed)
        error('IEle:InvalidInput','NodeSeed can not be a logical index')
    end
    
    if nargin< 5 || isempty(CreateFigure)
        CreateFigure=false;
    end
    
    if CtrlVar.TriNodes~=3
    
        error("EleFlooding:UnsupportedElementType","Currenlty only 3-node elements supported")
    
    end
    
    
    TRI=MUA.connectivity(EleSubset,:) ;
    
    % M = connectivity2adjacency(TRI)>0 ;
    % G=graph(M) ;
    
    G=graph(TRI,TRI(:,[2 3 1]));
    bins=conncomp(G) ;
    ID=bins(NodeSeed) ;
    nodes=find(bins==ID);
    IEle=any(ismember(MUA.connectivity,nodes),2);
    
    
    if CreateFigure
        figure ;
        hold off;
        PlotMuaMesh(CtrlVar,MUA,[],'k') ;
        hold on ; PlotMuaMesh(CtrlVar,MUA,IEle,'r') ;
    end
    
    
    
end
%%
