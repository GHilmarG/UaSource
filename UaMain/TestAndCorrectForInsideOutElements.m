function [connectivity,isChanged]=TestAndCorrectForInsideOutElements(CtrlVar,coordinates,connectivity,detJ)

%
% connectivity=TestAndCorrectForInsideOutElements(CtrlVar,coordinates,connectivity)
% Tests for inside-out elements and if any found, attempts to correct the connectivity
%
% If connectivity is changed, the logical variable isChanges is true on return. 

points=[1/3 1/3];
if nargin<4 || isempty(detJ)
    [~,detJ]=derivVector(coordinates,connectivity,1,points,1); % Nele x dof x nod
end

isChanged=1;

if all(detJ>0) 
    isChanged=0;
    return
end


% first check if all elements are inside-out as this quite a frequent situation and
% happens, for example, each time the
% orientation of the MeshBoundaryCoordinates is not consistent with gmsh/mesh2d

if all(detJ<0) 
    fprintf(' All elements are inside-out. All elements therefore flipped. \n  ') ;
    connectivity=FlipElements(connectivity);
    [~,detJ]=derivVector(coordinates,connectivity,1,points,1); % Nele x dof x nod
    if all(detJ<0) 
        fprintf(' All elements still inside-out. \n  ') ;
        save DumpFileTestAndCorrectForInsideOutElements
        error('Ua:TestAndCorrectForInsideOutElements:InsideOutelements','InsideOutElements')
    else
        fprintf(' No longer all elements inside-out. \n  ') ;
    end
    
end



if any(detJ<0)   %
    
    InsideOutElements=find(detJ<0);
    
    fprintf(' Out of %-i %-i-nod elements, %-i are inside-out.  ',size(connectivity,1),size(connectivity,2),numel(InsideOutElements)) ;
    fprintf('Attempting to correct for inside-out elements.')
    connectivity=FlipElements(connectivity,InsideOutElements);
    
    
    
    [~,detJ]=derivVector(coordinates,connectivity,1,points,1); % Nele x dof x nod
    if any(detJ<0) 
        InsideOutElements=find(detJ<0);
        fprintf('Out of %-i %-i-nod elements, %-i are still inside-out.\n  ',...
            size(connectivity,1),size(connectivity,2),numel(InsideOutElements)) ;
        save DumpFileTestAndCorrectForInsideOutElements
        error('Ua:TestAndCorrectForInsideOutElements:InsideOutelements','InsideOutElements')
    end
end


% are there still problems?

if all(detJ>0) 
    fprintf(' Successfully corrected for inside-out elements. \n  ') ;
else
    fprintf(' Despite best attempts, still some elements inside-out. \n  ') ;
    fprintf(' Nothing that can be done but to panic. No hope, give up, do not even try, just give up! \n  ') ;
    
    CtrlVar.PlotMesh=1;
    figure(1500) ; CtrlVar.MeshColor='k' ;
    PlotFEmesh(coordinates,connectivity,CtrlVar) ;
    hold on ;
    CtrlVar.MeshColor='r';
    CtrlVar.PlotLabels=1;
    CtrlVar.PlotNodes=1;
    
    
    PlotFEmesh(coordinates,connectivity,CtrlVar,InsideOutElements)
    title(' Inside-out elements in red')
    fprintf('[----  inside-out element connectivity: \n')
    for icounter=1:numel(InsideOutElements)
        fprintf('Element %i : ',InsideOutElements(icounter))
        fprintf('%i ', connectivity(InsideOutElements(icounter),:))
        fprintf('\n')
    end
    
    fprintf('-----------------------------------]    \n')
    save DumpFileTestAndCorrectForInsideOutElements
    error('Ua:TestAndCorrectForInsideOutElements:InsideOutelements','InsideOutElements')
    
end

end
