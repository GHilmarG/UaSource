function connectivity=TestAndCorrectForInsideOutElements(CtrlVar,coordinates,connectivity)

%
% connectivity=TestAndCorrectForInsideOutElements(CtrlVar,coordinates,connectivity)
% Tests for inside-out elements and if any found, attempts to correct the connectivity
%

[~,detJ]=derivVector(coordinates,connectivity,1,1); % Nele x dof x nod

% first check if all elements are inside out as this quite a frequent situation and
% happens, for example, each time the
% orientation of the MeshBoundaryCoordinates is not consistent with gmesh/mesh2d
if all(detJ<0) ;
    fprintf(' All elements are inside out. All elements therefore flipped. \n  ') ;
    connectivity=FlipElements(connectivity);
    [~,detJ]=derivVector(coordinates,connectivity,1,1); % Nele x dof x nod
    if all(detJ<0) ;
        fprintf(' All elements still inside out. \n  ') ;
        save TestSave
        error('TestAndCorrectForInsideOutElements: InsideOutElements')
    else
        fprintf(' No longer all elements inside out. \n  ') ;
    end
    
end


if any(detJ<0)   %
    
    InsideOutElements=find(detJ<0);
    
    fprintf(' Out of %-i %-i-nod elements, %-i are inside out.\n  ',size(connectivity,1),size(connectivity,2),numel(InsideOutElements)) ;
    fprintf('Correcting inside out elements. \n')
    connectivity=FlipElements(connectivity,InsideOutElements);
    
    
    [~,detJ]=derivVector(coordinates,connectivity,1,1); % Nele x dof x nod
    if any(detJ<0) ;
        InsideOutElements=find(detJ<0);
        fprintf('Out of %-i %-i-nod elements, %-i are still inside out.\n  ',...
            size(connectivity,1),size(connectivity,2),numel(InsideOutElements)) ;
        save TestSave
        error('TestAndCorrectForInsideOutElements: InsideOutElements')
    end
end


% are there still problems?

if all(detJ>0) ;
    fprintf(' No inside-out elements. \n  ') ;
else
    fprintf(' Despite best attempts, still some elements inside-out. \n  ') ;
    fprintf(' Nothing to do, but panic. No hope, give up, do not even try, just give up. \n  ') ;
    
    CtrlVar.PlotMesh=1;
    figure(1500) ; CtrlVar.MeshColor='k' ;
    PlotFEmesh(coordinates,connectivity,CtrlVar) ;
    hold on ;
    CtrlVar.MeshColor='r';
    CtrlVar.PlotLabels=1;
    CtrlVar.PlotNodes=1;
    
    
    PlotFEmesh(coordinates,connectivity,CtrlVar,InsideOutElements)
    title(' Inside out elements in red')
    fprintf('[----  inside-out element connectivity: \n')
    for icounter=1:numel(InsideOutElements)
        fprintf('Element %i : ',InsideOutElements(icounter))
        fprintf('%i ', connectivity(InsideOutElements(icounter),:))
        fprintf('\n')
    end
    
    
    fprintf('-----------------------------------]    \n')
    
    save TestSave
    error('TestAndCorrectForInsideOutElements: InsideOutElements')
end

end
