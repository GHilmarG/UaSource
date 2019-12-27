
function  I=MuaElementsContainingGivenNodes(CtlrVar,MUA,NodeList,EleList,AllOrAny)
    
    
    
    
    %%
    %   I=MuaElementsContainingGivenNodes(CtlrVar,MUA,NodeList,EleList,AllOrAny)
    %
    % NodeList must be an index variable, ie not a logical index
    %
    % Returns a logical list of elements containing either all of the nodes in
    % NodeList or one or more of the nodes in NodeList
    %
    % Optionally EleList can be specified, in which case the search is limited to the
    % elements in the list.  If EleList is not specified, the search is over all elements in
    % connectivity.
    %
    %  AllOrAny is a string variable and can be eitehr "all" or "any" depending
    %  on if elemetns containing all the nodes or any of the nodes should be
    %  found.
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
            
            if nargin<4  || isempty(EleList)
                I=any(ismember(MUA.connectivity,NodeList),2);
            else
                I=any(ismember(MUA.connectivity(EleList,:),NodeList),2);
            end
            
        case "all"
            
            if nargin<4 || isempty(EleList)
                I=all(ismember(MUA.connectivity,NodeList),2);
            else
                I=all(ismember(MUA.connectivity(EleList,:),NodeList),2);
            end
            
        otherwise
            error('sfad')
    end
    
end


