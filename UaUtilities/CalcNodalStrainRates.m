function [exxb,eyyb,exyb,exxd,eyyd,exyd]=CalcNodalStrainRates(CtrlVar,MUA,ub,vb,ud,vd)

%% Calculates strain rates at nodes
%
%
% Example
%
%    [exxb,eyyb,exyb,exxd,eyyd,exyd]=CalcNodalStrainRates(CtrlVar,MUA,ub,vb,ud,vd)
%
%   [exxb,eyyb,exyb]=CalcNodalStrainRates(CtrlVar,MUA,ub,vb)
%

if ~nargin==4 || ~nargin==6
    error('Ua:CalcNodalStrainRates','Wrong number of input arguments')
end

if ~isempty(ub)
    
    [dubdx,dubdy]=calcFEderivativesMUA(ub,MUA,CtrlVar);
    [dvbdx,dvbdy]=calcFEderivativesMUA(vb,MUA,CtrlVar);
    exxb=dubdx;
    eyyb=dvbdy;
    exyb=0.5*(dubdy+dvbdx);
    [exxb,eyyb,exyb]=ProjectFintOntoNodes(MUA,exxb,eyyb,exyb);
    
end

if nargin == 6
    if ~isempty(ud)
        [duddx,duddy]=calcFEderivativesMUA(ud,MUA,CtrlVar);
        [dvddx,dvddy]=calcFEderivativesMUA(vd,MUA,CtrlVar);
        
        exxd=duddx;
        eyyd=dvddy;
        exyd=0.5*(duddy+dvddx);
        
        [exxd,eyyd,exyd]=ProjectFintOntoNodes(MUA,exxd,eyyd,exyd);
    end
end


end


