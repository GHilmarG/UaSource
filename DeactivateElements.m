function  [coordinates,connectivity]=DeactivateElements(CtrlVar,ElementsToBeDeactivated,coordinates,connectivity)

% [coordinates,connectivity]=DeactivateElements(CtrlVar,iDeactivatedElements,coordinates,connectivity)
% Deactivates elements in the list iDeactivatedElements
% Nodes that are no longer part of the FE mesh are deleted and the connectivity updated accordingly
% Elements and nodes are renumbered

if ~any(ElementsToBeDeactivated)
    return
end

connectivity(ElementsToBeDeactivated,:)=[];

% eliminate coordinates that are no longer part of mesh, and update connectivity accordingly
[K,~,J]=unique(connectivity(:));
connectivity=reshape(J,size(connectivity));
coordinates=coordinates(K,:);

% renumber nodes and elements
if CtrlVar.sweep
    [coordinates,connectivity] = NodalSweep(coordinates,connectivity,CtrlVar.SweepAngle);
    [coordinates,connectivity] = ElementSweep(coordinates,connectivity,CtrlVar.SweepAngle);
end

end
