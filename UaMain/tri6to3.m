function [coordinates,connectivity]=tri6to3(coordinates,connectivity)
    
%%
%
%   [coordinates,connectivity]=tri6to3(coordinates,connectivity)
%
% creates 3 node triangles from 6 node tringles
    
    
    
    [Nele,nod]=size(connectivity);
    
    connectivity3=zeros(Nele,3);
    
    connectivity3(:,1)=connectivity(:,1);
    connectivity3(:,2)=connectivity(:,3);
    connectivity3(:,3)=connectivity(:,5);
    
    connectivity=connectivity3;
    
    [coordinates,connectivity]=RemoveNodesNotPartOfAnyElement(coordinates,connectivity);
    
    
    %CtrlVar.PlotNodes=1; CtrlVar.PlotLabels=1; figure ; PlotFEmesh(coordinates,connectivity,CtrlVar)
    
end

