
function [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)


%%
%
% User input m-file to define A and n in the Glenn-Steinemann flow law
%
%   [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)
%
% Usually A is defined on the nodes (but sometimes in an inverse run A might be
% defined as an element variable.)
% 
%% 

n=3 ;
A=6.338e-25;
AGlen=A*1e9*365.2422*24*60*60+zeros(MUA.Nnodes,1);


end

