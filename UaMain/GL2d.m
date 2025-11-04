function GF= GL2d(B,S,h,rhow,rho,connectivity,CtrlVar)

%% creates node and element masks with 1 if grounded, 0 otherwise
%
% GF = GL2d(B,S,h,rhow,rho,connectivity,CtrlVar)
%
% GF.node  : a vector with Nnodes values, 1 if node grounded, 0 otherwise
% GF.ele   : a vector with Nele values,   mean of nodal values
%
%

hf=(S-B)*rhow./rho ;  % hf is the maximum ice thickness possible for ice not to be grounded

GF.node= HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);
GF.ele=Nodes2EleMean(connectivity,GF.node);



end

