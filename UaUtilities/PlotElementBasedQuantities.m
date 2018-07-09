function [FigHandle,ColorbarHandle,tri]=PlotElementBasedQuantities(connectivity,coordinates,Value,CtrlVar,varargin)
    
% [FigHandle,ColorbarHandle,tri]=PlotElementBasedQuantities(connectivity,coordinates,Value,CtrlVar,varargin)
% plots elements scalar quantities given connectivity and coordinates
% does so using the matlab `patch' command
% CtrlVar is optional
% varargin is also optional and if given is simply passed on to the matlab patch command
% 
% Only good for plotting element quantities with one value per element
% (can not be used to plot element-based quantities defined at integration points)


    if size(coordinates,2)~=2
        error('coordinates in second argument')
    end
    
    [Nele,nod]=size(connectivity) ;
    Nvalues=length(Value);
    
        
    if Nele~=Nvalues
        error('Number of element-based values must be equal to number of elements.')
    end
    
    if nargin<4 || isempty(CtrlVar)
        CtrlVar.PlotXYscale=1;       
        CtrlVar.PlotsXaxisLabel=' ';
        CtrlVar.PlotsYaxisLabel=' ';
    end
    
    if nod==3
        tri=connectivity;
    elseif size(connectivity,2)==6
        tri=zeros(length(connectivity),3);
        tri(:,1)=connectivity(:,1);
        tri(:,2)=connectivity(:,3);
        tri(:,3)=connectivity(:,5);
    elseif size(connectivity,2)==10
        tri=zeros(length(connectivity),3);
        tri(:,1)=connectivity(:,1);
        tri(:,2)=connectivity(:,4);
        tri(:,3)=connectivity(:,7);
    end
    
    FigHandle=patch('Faces',tri,'Vertices',coordinates/CtrlVar.PlotXYscale,'FaceVertexCdata',Value,...
        'CDataMapping','scaled','FaceColor','flat','EdgeColor','none','FaceLighting','none',varargin{:}) ;
        
    ColorbarHandle=colorbar;
    xlabel(CtrlVar.PlotsXaxisLabel)  ; ylabel(CtrlVar.PlotsYaxisLabel) ;
    axis equal tight
    
    return
end