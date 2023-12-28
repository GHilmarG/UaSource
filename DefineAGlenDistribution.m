
function  [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,F)


%%
%
% User input m-file to define A and n in the Glenn-Steinemann flow law
%
%   [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,F)
% 
%   [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)
%
% Note: Use
%
%   [AGlen,B]=AGlenVersusTemp(T)
%
% to get A in the units kPa^{-3} yr^{-1} for some temperature T (degrees Celsius) 
%
%% 

n=3 ;
A=6.338e-25;                                           % A in SI units for temperate ice ; 
AGlen=A*1e9*365.2422*24*60*60+zeros(MUA.Nnodes,1);     % A in the units kPa^{-3} yr^{-1}


end

