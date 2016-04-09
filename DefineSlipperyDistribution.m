function [C,m]=DefineSlipperyDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

%%
%  User input m-file to define C and m in the Weertman sliding law
% [C,m]=DefineSlipperyDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)
%
% Usually C is defined on the nodes, but sometimes in an inverse run C might be
% defined as an element variable. The user makes this decision by setting
% CtrlVar.CisElementBased to true or false in Ua2D_InitialUserInput

warning('Ua:DefaultDefine','Using default DefineSlipperyDistribution \n')

m=3 ; C=0.01 ; 

if CtrlVar.CisElementBased
    C=C+zeros(MUA.Nele,1);
else
    C=C+zeros(MUA.Nnodes,1);
end

end
