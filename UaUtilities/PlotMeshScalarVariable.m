function [FigHandle,ColorbarHandle]=PlotMeshScalarVariable(CtrlVar,MUA,Variable,varargin)

% PlotMeshScalarVariable(CtrlVar,MUA,Variable,varargin)
% path plot of variable
% variable can be either nodal or element variable


persistent NodTri EleTri Nele Nnodes nod

[N,M]=size(Variable);

if N==MUA.Nnodes && M==1   % nodal variable
    
    if isempty(NodTri) || MUA.Nele~=Nele || MUA.Nnodes~= Nnodes || MUA.nod~=nod
        NodTri=MUA.connectivity;
    end
    
    [FigHandle,ColorbarHandle,NodTri]=PlotNodalBasedQuantities(NodTri,MUA.coordinates,Variable,CtrlVar,varargin{:});
    
elseif N==MUA.Nele && M==1 % element variable
    
    if isempty(EleTri) || MUA.Nele~=Nele || MUA.Nnodes~= Nnodes || MUA.nod~=nod
        EleTri=MUA.connectivity;
    end
    
    [FigHandle,ColorbarHandle,EleTri]=PlotElementBasedQuantities(EleTri,MUA.coordinates,Variable,CtrlVar,varargin{:});
    
else
    
    error('Ua:PlotMeshScalarVariable:InconsistentDimentionis','Inconsitent dimentions')
    
    
end

Nele=MUA.Nele ; Nnodes=MUA.Nnodes; nod=MUA.nod; 

end