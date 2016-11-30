

function [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

%%
%  User input m-file to define A and n in the Glenn-Steinemann flow law
%
%   [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)
%
% Usually A is defined on the nodes, but sometimes in an inverse run A might be
% defined as an element variable. The user makes this decision by setting
% CtrlVar.AGlenisElementBased to true or false in Ua2D_InitialUserInput
%%
warning('Ua:DefaultDefine','Using default DefineAGlenDistribution')


n=3 ;

if CtrlVar.AGlenisElementBased
    AGlen=3.0e-9+zeros(MUA.Nele,1) ;  % kPa year about -20 degrees Celcius
else
    AGlen=3.0e-9+zeros(MUA.Nnodes,1) ; % kPa year about -20 degrees Celcius
end


end

