function [FigHandle,ColorbarHandle]=PlotMeshScalarVariable(CtrlVar,MUA,Variable,varargin)

%%
% Plots a scalar variable over the mesh domain as a patch (see help patch)
%
% PlotMeshScalarVariable(CtrlVar,MUA,Variable,varargin)
%
% variable can be either nodal or element variable
%
%
% vararing is passed on to the patch command
%
% Examples:
%
% figure ; PlotMeshScalarVariable(CtrlVar,MUA,h)   
% plots the thickness distribution
%
%%

persistent NodTri EleTri Nele Nnodes nod

[N,M]=size(Variable);

if N==MUA.Nnodes && M==1   % nodal variable
    
    if isempty(NodTri)
        NodTri=MUA.connectivity;
    elseif MUA.Nele~=Nele || MUA.Nnodes~= Nnodes || MUA.nod~=nod 
        NodTri=MUA.connectivity;
    end
    
    [FigHandle,ColorbarHandle,NodTri]=PlotNodalBasedQuantities(NodTri,MUA.coordinates,Variable,CtrlVar,varargin{:});
    
elseif N==MUA.Nele && M==1 % element variable
    
    if isempty(EleTri) 
        EleTri=MUA.connectivity;
    elseif MUA.Nele~=Nele || MUA.Nnodes~= Nnodes || MUA.nod~=nod
        EleTri=MUA.connectivity;
    end
    
    [FigHandle,ColorbarHandle,EleTri]=PlotElementBasedQuantities(EleTri,MUA.coordinates,Variable,CtrlVar,varargin{:});
    
else
    
    fprintf('PlotMeshScalarVariable: Variable has inconsistent dimentions and can not be plotted.\n') 
    warning('Ua:PlotMeshScalarVariable:Inconsistentdimensions','Inconsistent dimensions')
    
    
end

Nele=MUA.Nele ; Nnodes=MUA.Nnodes; nod=MUA.nod; 

end