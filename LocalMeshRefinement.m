function [MUAnew,RunInfo]=LocalMeshRefinement(CtrlVar,RunInfo,MUAold,ElementsToBeRefined,ElementsToBeCoarsened)

persistent wasRefine



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


%%


narginchk(5,5)


if CtrlVar.InfoLevelAdaptiveMeshing>=1
    
    fprintf('\t Local adaptive remeshing at run-step %i and time=%f \n',CtrlVar.CurrentRunStepNumber,CtrlVar.time)
    
end





MUAnew=MUAold;
MUAonInput=MUAold;
% Make sure that lists are logical
if ~islogical(ElementsToBeRefined)
    
    T= false(MUAold.Nele,1);
    T(ElementsToBeRefined)=true;
    ElementsToBeRefined=T;
    
end

if ~islogical(ElementsToBeCoarsened)
    
    T= false(MUAold.Nele,1);
    T(ElementsToBeCoarsened)=true;
    ElementsToBeCoarsened=T;
    
end



nRefine=numel(find(ElementsToBeRefined));
nCoarsen=numel(find(ElementsToBeCoarsened));




if nRefine <CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan  && nCoarsen < CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan
    
    MUAnew=MUAold;
    
    if CtrlVar.InfoLevelAdaptiveMeshing>=1
        if nRefine==0 && nCoarsen==0
            fprintf('LocalMeshRefinement: No elements to be refined or coarsened.\n')
        else
            fprintf('LocalMeshRefinement: The numbers of elements to be refined (%i) or coarsened (%i) are both less than CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan (%i) \n',...
                nRefine,nCoarsen,CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan)
            fprintf('                     Hence, no local mesh refinement is done.\n')
        end
    end
    
    return
end

if isempty(wasRefine)
    wasRefine=0;
end


% refine and smoothmesh only works for 3-nod elements
% so first change to 3-nod if needed
[MUAold.coordinates,MUAold.connectivity]=ChangeElementType(MUAold.coordinates,MUAold.connectivity,3);

switch CtrlVar.MeshRefinementMethod
    
    case {'explicit:local:red-green'}
        
        
        RunInfo.MeshAdapt.Method='Red-Green Refinement';
        [MUAold.coordinates,MUAold.connectivity] = refine(MUAold.coordinates,MUAold.connectivity,ElementsToBeRefined);
        
        if CtrlVar.LocalAdaptMeshSmoothingIterations>0
            opts=CtrlVar.Smooth2.opts;
            tnum=[]; edge=[];
            [MUAold.coordinates,edge,MUAold.connectivity,tnum] = smooth2(MUAold.coordinates,edge,MUAold.connectivity,tnum,opts);
        end
        %        [MUAold.coordinates] = GHGsmoothmesh(MUAold.coordinates,MUAold.connectivity,CtrlVar.LocalAdaptMeshSmoothingIterations,[]);
        %        MUAold.connectivity=FlipElements(MUAold.connectivity);
        
        % create MUA (this takes care of any eventual change in element type as well)
        MUAnew=CreateMUA(CtrlVar,MUAold.connectivity,MUAold.coordinates);
        
        
    case 'explicit:local:newest vertex bisection'
        
        %MUAold.connectivity=FlipElements(MUAold.connectivity);
        if  ~isfield(MUAold,'RefineMesh')  ||  isempty(MUAold.RefineMesh)
            mesh = genMesh(MUAold.connectivity, MUAold.coordinates);
            mesh.bd=[];
            mesh = genBisectionMesh(mesh);
            mesh = SelectRefinementEdge(mesh);
            MUAold.RefineMesh=mesh;
        else
            mesh=MUAold.RefineMesh;
        end
        
        
        if MUAold.nod~=3
            mesh.TR=triangulation(mesh.elements,mesh.coordinates);
        end
        
        
        Na=size(mesh.coordinates,1);  Ea=size(mesh.elements,1);
        
        % do refinement/coarsening alternativily, unless if there are
        % significant in number of elements to be refined/coarsened
        
        
        if wasRefine  &&  ( nRefine < 2*nCoarsen) || nRefine==0
            isRefine=0;
        else
            isRefine=1;
        end
        
        
        wasRefine=isRefine;
        
        if isRefine
            RunInfo.MeshAdapt.Method='Bisection Refinement';
            fprintf(' Refining %i elements \n',nRefine)
            if MUAold.nod~=3
                meshElementsToBeRefined=pointLocation(mesh.TR,[MUAold.xEle(ElementsToBeRefined) MUAold.yEle(ElementsToBeRefined)]);
            else
                meshElementsToBeRefined=ElementsToBeRefined;
            end
            
            mesh = bisectionRefine2D(mesh,meshElementsToBeRefined);
            
        else
            RunInfo.MeshAdapt.Method='Bisection Coarsening';
            fprintf(' Coarsening %i elements \n',nCoarsen)
            %mesh.elements=FlipElements(mesh.elements);
            if MUAold.nod~=3
                meshElementsToBeCoarsened=pointLocation(mesh.TR,[MUAold.xEle(ElementsToBeCoarsened) MUAold.yEle(ElementsToBeCoarsened)]);
            else
                meshElementsToBeCoarsened=ElementsToBeCoarsened;
            end
            mesh = bisectionCoarsen(mesh,meshElementsToBeCoarsened);
        end
        
        %elements=TestAndCorrectForInsideOutElements(CtrlVar,mesh.coordinates,mesh.elements);
        
        Nb=size(mesh.coordinates,1);  Eb=size(mesh.elements,1);
        isMeshChangedSufficiently=abs(Eb-Ea) > CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan ;
        
        
        
        if isMeshChangedSufficiently
            fprintf('In local mesh-refinement step the change in the number of elements and nodes was %i and %i, respectivily  (#R/#C)=(%i/%i). \n',...
                Eb-Ea,Nb-Na,nRefine,nCoarsen)
            % This will also change the element type back to 6 adn 10 nodes if required
            MUAnew=CreateMUA(CtrlVar,mesh.elements,mesh.coordinates,mesh);
            
        else
            
            if CtrlVar.InfoLevelAdaptiveMeshing>=10
                
                if isRefine
                    fprintf('LocalMeshRefinement: The change in numbers of elements during refinement (%i) is less than CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan (%i) \n',...
                        Eb-Ea,CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan)
                    fprintf('                     Hence, this local mesh refinement is discarded and the mesh is not updated.\n')
                else
                    fprintf('LocalMeshRefinement: The change in numbers of elements during coarsening (%i) is less than CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan (%i) \n',...
                        Eb-Ea,CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan)
                    fprintf('                     Hence, this local mesh coarsening is discarded and the mesh is not updated.\n')
                end
                
            end
            
            fprintf('Mesh not changed in local mesh-refinement step (#R/#C)=(%i/%i). \n',nRefine,nCoarsen)
            
            if CtrlVar.ManuallyDeactivateElements && ~(size(mesh.elements,1)==MUAnew.Nele && size(mesh.coordinates,1)==MUAnew.Nnodes)
                % If the user wants to manually deactivate elements, I must
                % re-introduce the full mesh at each mesh refinement stage to
                % allow for the re-activation of regions previously deactivated
                %
                fprintf('LocalMeshRefinement: Because CtrlVar.ManuallyDeactivateElements is set to true, the initial mesh is now reintroduced. \n')
                fprintf('                     MUA is now again the mesh at start of run. \n')
                MUAnew=CreateMUA(CtrlVar,mesh.elements,mesh.coordinates,mesh);
            else
                MUAnew=MUAonInput;
            end
            
        end
        %%
        
        
    otherwise
        
        
        error('case not found.')
        
end


end