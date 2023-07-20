
function  I=AllElementsContainingGivenNodes(connectivity,NodeList,EleList)

    
    
    
%%
%   I=AllElementsContainingGivenNodes(connectivity,NodeList,EleList)
%
% NodeList must be an index variable, ie not a logical index
%
% Returns a logical list of elements containing one or more of the nodes in NodeList
%
% Optionally EleList can be specified, in which case the search is limited to the
% elements in the list.  If EleList is not specified, the search is over all elements in
% connectivity.
%
% Example: find all elements containing one or more of the  nodes 1, 10 , 11
%
%
%  I=AllElementsContainingGivenNodes(MUA.connectivity,[1 10 11]); 
%  figure ; PlotMuaMesh(CtrlVar,MUA,[],'k') ;  hold on ; PlotMuaMesh(CtrlVar,MUA,I,'r') ;
%
%
%
%%

if islogical(NodeList)
    NodeList=find(NodeList) ;
end

if nargin<3
    I=any(ismember(connectivity,NodeList),2);
else
    I=any(ismember(connectivity(EleList,:),NodeList),2);
end

end


