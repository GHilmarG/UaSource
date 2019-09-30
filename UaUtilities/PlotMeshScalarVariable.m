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

persistent NodTri EleTri Nele Nnodes nod

if islogical(Variable)
    
    Variable=double(Variable);
    
end


[N,M]=size(Variable);

if N==MUA.Nnodes && M==1   % nodal variable
    
    if isempty(NodTri) || isempty(Nnodes) 
        NodTri=MUA.connectivity;
    elseif MUA.Nele~=Nele || MUA.Nnodes~= Nnodes || MUA.nod~=nod 
        NodTri=MUA.connectivity;
    end
    
    [FigHandle,ColorbarHandle,NodTri]=PlotNodalBasedQuantities(NodTri,MUA.coordinates,Variable,CtrlVar,varargin{:});
    
elseif N==MUA.Nele && M==1 % element variable
    
    if isempty(EleTri) || isempty(Nele) 
        EleTri=MUA.connectivity;
    elseif MUA.Nele~=Nele || MUA.Nnodes~= Nnodes || MUA.nod~=nod
        EleTri=MUA.connectivity;
    end
    
    [FigHandle,ColorbarHandle,EleTri]=PlotElementBasedQuantities(EleTri,MUA.coordinates,Variable,CtrlVar,varargin{:});
    
else
    
    fprintf('PlotMeshScalarVariable: Variable has inconsistent dimensions and can not be plotted.\n') 
    warning('Ua:PlotMeshScalarVariable:Inconsistentdimensions','Inconsistent dimensions')
    
    
end

Nele=MUA.Nele ; Nnodes=MUA.Nnodes; nod=MUA.nod; 

end