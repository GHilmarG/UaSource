function [coordinates,connectivity]=tri10to3(coordinates,connectivity)
    
    % creates 3 node triangles from 10 node triangle
    
    
    
    [Nele,nod]=size(connectivity);
    
    connectivity3=zeros(Nele,3);
    
    connectivity3(:,1)=connectivity(:,1);
    connectivity3(:,2)=connectivity(:,4);
    connectivity3(:,3)=connectivity(:,7);
    
    
    connectivity=connectivity3;
    
    [coordinates,connectivity]=RemoveNodesNotPartOfAnyElement(coordinates,connectivity);
    
    
    %CtrlVar.PlotNodes=1; CtrlVar.PlotLabels=1; figure ; PlotFEmesh(coordinates,connectivity,CtrlVar)
    
end