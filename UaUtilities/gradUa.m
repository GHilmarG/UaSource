function [dfdx,dfdy]=gradUa(CtrlVar,MUA,f)


[dfdx,dfdy]=calcFEderivativesMUA(f,MUA,CtrlVar) ; 
[dfdx,dfdy]=ProjectFintOntoNodes(CtrlVar,MUA,dfdx,dfdy) ;

end
