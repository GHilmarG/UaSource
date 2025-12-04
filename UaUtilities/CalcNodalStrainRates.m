function [exxb,eyyb,exyb,exxd,eyyd,exyd]=CalcNodalStrainRates(MUA,ub,vb,ud,vd)


%%
%  
% Calculates horizontal strain rates given horizontal velocities. 
%
%    [exxb,eyyb,exyb,exxd,eyyd,exyd]=CalcNodalStrainRates(MUA,ub,vb,ud,vd)
%
% Returns nodal values.
%
% Note: The strain rates are calculated at integration points and then projected
% onto nodes. 
%
% The projection does not conserve positivity and positve integration
% values can become negative at nodes. The effectiv strain rate , e, is for
% this reason calculated directly from nodal values, ensuring that e is
% always positive.
%
% Example: 
%
% Read data, calculate nodal strain rates, and then plot over FE mesh at roughly equal spaced grid.
%
%   load ('GaussPeak_Example_Restartfile.mat','MUA','CtrlVarInRestartFile','F','GF','BCs');  % load data
%   CtrlVar=CtrlVarInRestartFile; x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2);
%   [exx,eyy,exy]=CalcNodalStrainRates(CtrlVar,MUA,F.ub,F.vb);                             % calculate strain rates
%   [X,Y]=ndgrid(linspace(min(x),max(x),10),linspace(min(y),max(y),10));
%   I=nearestNeighbor(MUA.TR,[X(:) Y(:)]);  % find nodes within computational grid closest to the regularly scape X and Y grid points.
%   figure
%   CtrlVar.PlotNodes=0; PlotMuaMesh(CtrlVar,MUA,[],'color','k') ;                           % Plot FE mesh
%   hold on
%   scale=1e5; LineWidth=2 ;
%   PlotTensor(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,exx(I),exy(I),eyy(I),scale,LineWidth);  % plot strain rates
%   axis equal
%
%
% See also CalcNodalStrainRatesAndStresses

narginchk(3,5)

if nargin==4 
       fprintf("The parameter list no longer contains CtrlVar \n")
end

if ~nargin==3 || ~nargin==5
    error('Ua:CalcNodalStrainRates','Wrong number of input arguments')
end

if ~isempty(ub)
    
    [dubdx,dubdy]=calcFEderivativesMUA(ub,MUA);
    [dvbdx,dvbdy]=calcFEderivativesMUA(vb,MUA);
    exxb=dubdx;
    eyyb=dvbdy;
    exyb=0.5*(dubdy+dvbdx);
    [exxb,eyyb,exyb]=ProjectFintOntoNodes(MUA,exxb,eyyb,exyb);
else
    exxb=[] ; eyyb=[] ; exyb=[]; 
end

if nargin == 5
    if ~isempty(ud)
        [duddx,duddy]=calcFEderivativesMUA(ud,MUA);
        [dvddx,dvddy]=calcFEderivativesMUA(vd,MUA);
        
        exxd=duddx;
        eyyd=dvddy;
        exyd=0.5*(duddy+dvddx);
        
        [exxd,eyyd,exyd]=ProjectFintOntoNodes(MUA,exxd,eyyd,exyd);
    end
else
    exxd=[] ; eyyd=[] ; exyd=[]; 
end


end


