function [FigHandle,ColorbarHandle]=PlotMeshScalarVariable(CtrlVar,MUA,Variable,varargin)

%%
% Plots a scalar variable over the mesh domain as a patch (see help patch)
%
% PlotMeshScalarVariable(CtrlVar,MUA,Variable,varargin)
%
% variable can be either nodal, element variable (i.e. single value per element) or an integration-point variable.
%
% If variable is empty, nothing is plotted and no warning given.
%
% If 'Variable' not empty, but is not nodal, element or integration-point variable, I complain a bit.
%
% varargin is passed on to the patch command
%
% *Examples:*
%
% Plot element sizes
%
%   load('MUA-PIG-TWG-Example.mat','MUA','BCs','CtrlVar')
%   Tarea=TriAreaFE(MUA.coordinates,MUA.connectivity); Tlength=sqrt(2*Tarea) ;
%   figure ; PlotMeshScalarVariable(CtrlVar,MUA,Tlength) ; title('Element sizes')
%
%
% Plot a nodal variable (here as an example, the x coordinates of the nodes)
%
%   load('MUA-PIG-TWG-Example.mat','MUA','BCs','CtrlVar')
%   x=MUA.coordinates(:,1);
%   figure ; PlotMeshScalarVariable([],MUA,x) ;
%
% Plot the floating mask:
%
%   load('MUA-PIG-TWG-Example.mat','MUA','GF','CtrlVar')
%   x=MUA.coordinates(:,1);
%   figure ; PlotMeshScalarVariable(CtrlVar,MUA,GF.node) ; title('The nodal floating mask (floating=0, grounded=1)')
%%



if islogical(Variable)

    Variable=double(Variable);

end


[N,M]=size(Variable);
FigHandle=[] ;
ColorbarHandle=[] ;

if N==MUA.Nnodes && M==1   % nodal variable

    NodTri=MUA.connectivity;

    [FigHandle,ColorbarHandle,NodTri]=PlotNodalBasedQuantities(NodTri,MUA.coordinates,Variable,CtrlVar,varargin{:});

elseif N==MUA.Nele && M==1 % element variable


    EleTri=MUA.connectivity;
    [FigHandle,ColorbarHandle,EleTri]=PlotElementBasedQuantities(EleTri,MUA.coordinates,Variable,CtrlVar,varargin{:});

elseif N==MUA.Nele && M==MUA.nip % integration-point  variable


    % This case is slighly more complicated, because the set of integration point can have duplicates if integration points fall on the
    % element edges, and one must also get rid of any resulting triangles outside of (a possible non-convex) domain.



    x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
    [xint,yint] = CalcIntegrationPointsCoordinates(MUA);

    % create vectors Xint and Yint of unique integration points and triangulise that set of points
    Xint=xint(:) ; Yint=yint(:); [~, Iint, ~] = unique([Xint Yint],'first','rows'); Iint = sort(Iint); Xint = Xint(Iint); Yint = Yint(Iint);
    DTint = delaunayTriangulation(Xint,Yint);

    % get rid of triangles outside of the polygon define by MeshBoundaryCoordinates
    ic=incenter(DTint);
    [cnInt,on] = inpoly2(ic,[x(MUA.Boundary.EdgeCornerNodes) y(MUA.Boundary.EdgeCornerNodes)]);
    DTintTriInside=DTint.ConnectivityList(cnInt,:);



    [FigHandle,ColorbarHandle]=PlotIntegrationPointBasedQuantities(CtrlVar,DTintTriInside,DTint.Points,Variable,varargin{:}) ;



elseif ~isempty(Variable)

    fprintf('PlotMeshScalarVariable: Variable has inconsistent dimensions and can not be plotted.\n')
    warning('Ua:PlotMeshScalarVariable:Inconsistentdimensions','Inconsistent dimensions')
    FigHandle=[] ;
    ColorbarHandle=[] ;

end



end