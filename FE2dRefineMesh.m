
function [coordinates,connectivity]=FE2dRefineMesh(coordinates,connectivity)
    
    % refines mesh by subdividing all triangles into four
    
    [NeleIn,nodIn]=size(connectivity);
    
    
    % change all elements into 3-node elemetns
    [coordinates,connectivity]=ChangeElementType(coordinates,connectivity,3);
    
    % refine all 3-node elements by globally splitting into four elements
    [coordinates,connectivity] = refine(coordinates,connectivity);
    connectivity=FlipElements(connectivity);  % `refine.m' assumes different orientation of element, so I must change it here
    
    % change back to original element type
    [coordinates,connectivity]=ChangeElementType(coordinates,connectivity,nodIn);
    
    % sweep and renumber
    alphaSweep=0.01;
    [coordinates,connectivity] = ElementSweep(coordinates,connectivity,alphaSweep);
    [coordinates,connectivity] = NodalSweep(coordinates,connectivity,alphaSweep);
    
    
end
