function MUAnew=LocalMeshRefinement(CtrlVar,MUAold,ElementsToBeRefined)

% Local mesh refinement by subdividing selected triangular elements into four elements
% MUAnew=LocalMeshRefinement(CtrlVar,MUAold,ElementsToBeRefined)
%
% ElementsToBeRefined : either a logical vector with same number of elements as there are elements in the mesh,
%                       that is numel(ElementsToBeRefined)=MUA.Nele,
%                       or a list of numbers.
%
% Example:  To refine elements 200 and 500 in an existing mesh
%           MUAnew=LocalMeshRefinement(CtrlVar,MUAold,[200 500])
%
%
%          To refine all elements within a given radius R
%          xEle=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,1));
%          yEle=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,2));
%          r=sqrt(xEle.*xEle+yEle.*yEle);
%          I=r<R;
%          MUAnew=LocalMeshRefinement(CtrlVar,MUA,I)


% refine and smoothmesh only works for 3-nod elements
% so first change to 3-nod if needed

if ~islogical(ElementsToBeRefined)
    
    T= false(MUAold.Nele,1);
    T(ElementsToBeRefined)=true;
    
else
    
    T=ElementsToBeRefined;
    
end


[MUAold.coordinates,MUAold.connectivity]=ChangeElementType(MUAold.coordinates,MUAold.connectivity,3);
[MUAold.coordinates,MUAold.connectivity] = refine(MUAold.coordinates,MUAold.connectivity,T);
[MUAold.coordinates] = GHGsmoothmesh(MUAold.coordinates,MUAold.connectivity,CtrlVar.LocalAdaptMeshSmoothingIterations,[]);
MUAold.connectivity=FlipElements(MUAold.connectivity);

% create MUA (this takes care of any eventual change in element type as well)
MUAnew=CreateMUA(CtrlVar,MUAold.connectivity,MUAold.coordinates);



end