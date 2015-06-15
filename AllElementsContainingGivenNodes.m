
function I=AllElementsContainingGivenNodes(connectivity,NodeList,EleList)

%I=AllElementsContainingGivenNodes(connectivity,NodeList,EleList)
%
% Returns a logical list of elements containing one or more of the nodes in NodeList
%
% Optionally EleList can be specified, in which case the search is limited to the
% elements in the list.  If EleList is not specified, the search is over all elements in
% connectivity.
%
% Example: find all elements containing one or more of the  nodes 1, 10 , 11
%  I=AllElementsContainingGivenNodes(MUA.connectivity,[1 10 11]); 
%  find(I)                % lists the element numbers
%  MUA.connectivity(I,:)  % gives the connectivity of those elements
%  figure ; PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar,I) % plots the elements

if nargin<3
    I=any(ismember(connectivity,NodeList),2);
else
    I=any(ismember(connectivity(EleList,:),NodeList),2);
end

end


