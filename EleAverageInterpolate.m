function AveragedNodalValues=EleAverageInterpolate(NodalValues,coordinates,connectivity,TriAreas,CtrlVar)
    
    %% Averages from nodes to elements and then maps back to nodes
    %
    % AveragedNodalValues=EleAverageInterpolate(NodalValues,coordinates,connectivity,TriAreas,CtrlVar)
    %  
    %  TirAreas and CtrlVar are optional
    %
    %
    %
    
    persistent DT
    
    % Averages nodal values over each element
    % assigns that value to the center point of each triangle
    % and then interpolates back to nodes
    
    if nargin==3 
        CtrlVar.NormalizeWithAreas=0;
    end
    
    [Nele,nod]=size(connectivity);
    EleAverage=real(mean(reshape(NodalValues(connectivity),Nele,nod),2));
    xEle=mean(reshape(coordinates(connectivity,1),Nele,nod),2);
    yEle=mean(reshape(coordinates(connectivity,2),Nele,nod),2);
    
    % normalize with area?
    if CtrlVar.NormalizeWithAreas
        fprintf(' Normalizing ele averages with ele areas \n ')
        EleAverage=EleAverage./TriAreas;
    end
    
    
    if isempty(DT) || length(DT.X)~=length(xEle)
        fprintf('DelaunayTri in AveragedNodalValues \n ')
        DT = DelaunayTri(xEle,yEle) ;
    end
    
    
    AveragedNodalValues=Grid1toGrid2(DT,EleAverage,coordinates(:,1),coordinates(:,2));
    
 
    
    
end