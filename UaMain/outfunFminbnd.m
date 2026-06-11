function [stop,Outs]=outfunFminbnd(x,optimValues,state)

persistent xVector fVector iCount

if nargout==2 && nargin==0 
    Outs.x=xVector;
    Outs.f=fVector;
    stop=false ;
    iCount=[] ; % reset
    return
end


if isempty(iCount)  || optimValues.funccount==1 
    iCount=0 ;
    xVector=nan(100,1);
    fVector=nan(100,1);
end

stop=false;

iCount=iCount+1;

if ~isempty(x) && numel(x)==1 
    xVector(iCount) = x ;
    fVector(iCount) = optimValues.fval ;
end




end