

function [P,N]=CalcLevelSetPertubationFunctional(CtrlVar,MUA,f)


txt=char(CtrlVar.LevelSetFABCostFunction);

idp = strfind(CtrlVar.LevelSetFABCostFunction,"p");
idq = strfind(CtrlVar.LevelSetFABCostFunction,"q");

p=str2double(txt(idp+1:idp+1));
q=str2double(txt(idq+1:idq+1));

xa=CtrlVar.LSFslope;

[dfdx,dfdy]=calcFEderivativesMUA(f,MUA,CtrlVar);

N=sqrt(dfdx.*dfdx+dfdy.*dfdy);

fint=(N.^q-xa).^p;

Int=FEintegrate2Dint(CtrlVar,MUA,fint);
P=sum(sum(Int))/(p*q);

if nargout>1
    N=ProjectFintOntoNodes(MUA,N) ;
end

end


