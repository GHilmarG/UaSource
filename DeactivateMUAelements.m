function  MUA=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated)

%
% Deactivates elements in the list iDeactivatedElements
% Nodes that are no longer part of the FE mesh are deleted and the connectivity updated accordingly
% Elements and nodes are renumbered

if ~any(ElementsToBeDeactivated)
    return
end

MUA.connectivity(ElementsToBeDeactivated,:)=[];

% eliminate MUA.coordinates that are no longer part of mesh, and update MUA.connectivity accordingly
[K,~,J]=unique(MUA.connectivity(:));
MUA.connectivity=reshape(J,size(MUA.connectivity));
MUA.coordinates=MUA.coordinates(K,:);

% renumber nodes and elements
if CtrlVar.sweep
    [MUA.coordinates,MUA.connectivity] = NodalSweep(MUA.coordinates,MUA.connectivity,CtrlVar.SweepAngle);
    [MUA.coordinates,MUA.connectivity] = ElementSweep(MUA.coordinates,MUA.connectivity,CtrlVar.SweepAngle);
end

MUA=UpdateMUA(CtrlVar,MUA);

% MUA.RefineMesh=[] ; %  As I now have a new mesh I need to reset the newest vertex bisection data structure.
%                    %  Therefore no further unrefinement over the previous mesh can be done.



end
