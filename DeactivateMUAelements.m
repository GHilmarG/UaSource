function  [MUA,K]=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated)

%
% Deactivates elements in the list iDeactivatedElements
% Nodes that are no longer part of the FE mesh are deleted and the connectivity updated accordingly
% Elements and nodes are renumbered

if ~any(ElementsToBeDeactivated)
    K=[];
    return
end

MUA.connectivity(ElementsToBeDeactivated,:)=[];

% eliminate MUA.coordinates that are no longer part of mesh, and update MUA.connectivity accordingly
[K,~,ic]=unique(MUA.connectivity(:));
MUA.connectivity=reshape(ic,size(MUA.connectivity));
MUA.coordinates=MUA.coordinates(K,:);

% K is the subset of nodes that I keep.
%
% Assuming there are no further changes to the nodal numbering, if I wanted to interpolate, I could do xNew=xOld(k) ;

% renumber nodes and elements
if CtrlVar.sweep
    [MUA.coordinates,MUA.connectivity,p] = NodalSweep(MUA.coordinates,MUA.connectivity,CtrlVar.SweepAngle);
    [MUA.coordinates,MUA.connectivity] = ElementSweep(MUA.coordinates,MUA.connectivity,CtrlVar.SweepAngle);
    K=K(p) ; 
end

if CtrlVar.UpdateMUAafterDeactivating
    MUA=UpdateMUA(CtrlVar,MUA);
end


% MUA.RefineMesh=[] ; %  As I now have a new mesh I need to reset the newest vertex bisection data structure.
%                    %  Therefore no further unrefinement over the previous mesh can be done.



end
