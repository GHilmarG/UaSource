function  [MUA,k,l]=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated,kIn,lIn)

%%
%
% Deactivates elements in the list iDeactivatedElements
%
% The variable ElementsToBeDeactivated can be either a logical or an index array.
%
%
% Nodes that are no longer part of the FE mesh are deleted and the connectivity updated accordingly
%
% 
% If CtrlVar.sweep is true, elements and nodes are then renumbered.
%
%  k gives me the mapping between new and old node numbers, that is:
% 
%       j=k(i) 
% 
% gives the old node number j for the new node number i
%
% If needed, one can transfer/map the old nodal values, fOld, onto the new subset of nodes using:
%
%     fNew=fOld(k) ;
%
% and
%
%   fNew(i)=fOld(k(i)) ;
%
% where fNew are the nodal values over the new mesh MUA (ie MUA returned).
%
%  l gives the mapping between old and new nodes numbers, that is:
%
%        i=l(j) 
%
% gives the new node number i for the old node number j.
%
%   l(j) is nan where there is no node in the new mesh corresponding to the old node j, ie for all deleted nodes.
%
% To get a list of deleted old node numbers use:
%
%   find(isnan(l)) 
%
% And to get a logical list of the nodes on the old mesh that were kept in the new:
%
%   ~isnan(l)             ;   % logical list of nodes in the old mesh used/kept in the new mesh.
%
% If elements are deactivated repeatedly and one needs to know the mapping between the original and the final mesh do:
%
%
%   k=k1(k2(k3))     % where k1 results from the first, k2 from the second, etc, deactivations
%
% and then:
%
%     l=1:nNodesIn ; l=l(:)+nan;      % l gives the mapping betwenen old and new, i=l(j) gives the new node number i for the old node number j
%     l(k(1:numel(k)))=1:numel(k);    %  i=l(j) is nan if the original j node was deleted
%
% where nNodesIn is the number of nodes in the initial mesh, i.e. before the first round of deactivations.
%
%%

% This works equally for both logical and index arrays. 
if ~any(ElementsToBeDeactivated)
    k=1:MUA.Nnodes; 
    l=1:MUA.Nnodes; 
    return
end

nNodesIn=MUA.Nnodes;




MUA.connectivity(ElementsToBeDeactivated,:)=[];

% eliminate MUA.coordinates that are no longer part of mesh, and update MUA.connectivity accordingly
[k,~,ic]=unique(MUA.connectivity(:));
MUA.connectivity=reshape(ic,size(MUA.connectivity));
MUA.coordinates=MUA.coordinates(k,:);


% K is the subset of nodes that I keep.
%
% Assuming there are no further changes to the nodal numbering, if I wanted to interpolate, I could do xNew=xOld(k) ;

% renumber nodes and elements
if CtrlVar.sweep
    [MUA.coordinates,MUA.connectivity,p] = NodalSweep(MUA.coordinates,MUA.connectivity,CtrlVar.SweepAngle);
    [MUA.coordinates,MUA.connectivity] = ElementSweep(MUA.coordinates,MUA.connectivity,CtrlVar.SweepAngle);
    k=k(p) ;
end

if CtrlVar.UpdateMUAafterDeactivating
    MUA=UpdateMUA(CtrlVar,MUA);
else
    MUA.M=[];
    MUA.Deriv=[];
    MUA.DetJ=[];
    MUA.TR=[];
    MUA.Boundary=[];
    MUA.dM=[];
    MUA.xEle=[];
    MUA.yEla=[];
    MUA.Nnodes=size(MUA.coordinates,1);
    MUA.Nele=size(MUA.connectivity,1);
    MUA.EleAreas=TriAreaFE(MUA.coordinates,MUA.connectivity); % areas for each element
    MUA.Area=sum(MUA.EleAreas);

end

if nargin>4  && ~isempty(kIn)
    k=kIn(k);
end



if nargout==3   % create a mapping from old to new node numbers
    if nargin==5  && ~isempty(lIn)
        l=lIn+nan;
    else
        %l=1:nNodesIn ; l=l(:)+nan;
        l=nan(nNodesIn,1) ;
    end
    l(k(1:numel(k)))=1:numel(k);
end


% MUA.RefineMesh=[] ; %  As I now have a new mesh I need to reset the newest vertex bisection data structure.
%                    %  Therefore no further unrefinement over the previous mesh can be done.



end
