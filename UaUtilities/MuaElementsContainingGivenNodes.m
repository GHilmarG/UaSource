
function  I=MuaElementsContainingGivenNodes(CtrlVar,MUA,NodeList,EleList,AllOrAny)
    
    
    
    
    %%
    %   I=MuaElementsContainingGivenNodes(CtrlVar,MUA,NodeList,EleList,AllOrAny)
    %
    % Finds elements in *MUA* containing nodes in *NodeList* .
    %
    %
    % *NodeList* must be an index variable (indexing with element position), ie not a logical index
    %
    %
    % On return *I* is a logical index list of elements containing either ALL of the nodes in
    % *NodeList* or ANY (default) of the nodes in *NodeList*
    %
    % Optionally *EleList* can be specified, in which case the search is limited to the
    % elements in that list.
    %
    % *EleList* must be a logical index 
    %
    % *AllOrAny* is a string variable and can be either "all" or "any" (default) depending
    % on if elements containing all the nodes or any of the nodes should be
    % found.
    %
    % Example: find elements containing one or more nodes with thickness
    % less that 10, and plot those over the mesh.
    %
    %   NodeListLogical=  F.h <= 10;
    %   IEle=MuaElementsContainingGivenNodes(CtrlVar,MUA,find(NodeListLogical)) ;
    %   figure ; PlotMuaMesh(CtrlVar,MUA,[],'k') ;  hold on ; PlotMuaMesh(CtrlVar,MUA,IEle,'r') ;
    %
    % Example: find elements, within the element numbers ranging from 10 to
    % 100, where all nodes have a thickness of 10 or less.
    %
    %   NodeListLogical=  F.h <= 10;
    %   IEle=MuaElementsContainingGivenNodes(CtrlVar,MUA,find(NodeListLogical),10:100,"all") ;
    %   figure ; PlotMuaMesh(CtrlVar,MUA,[],'k') ;  hold on ; PlotMuaMesh(CtrlVar,MUA,IEle,'r') ;
    %
    %
    %%
    
    narginchk(3,5)
    
    if nargin<5 || isempty(AllOrAny)
        AllOrAny="any" ;
    end
    
    
    switch lower(AllOrAny)
        
        case "any"
            I=any(ismember(MUA.connectivity,NodeList),2);
            
        case "all"
            
            I=all(ismember(MUA.connectivity,NodeList),2);
            
        otherwise
            error('sfad')
    end
    
    if nargin >=4 && ~isempty(EleList)
        I=I & EleList ;
    end
    
end


