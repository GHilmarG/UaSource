function [UserVar,rho,rhow,g]=DefineDensities(UserVar,CtrlVar,MUA,time,s,b,h,S,B)

%%
% Define ice and ocean densities, as well as the gravitational acceleration.
%
%   [UserVar,rho,rhow,g]=DefineDensities(UserVar,CtrlVar,MUA,time,s,b,h,S,B)
%
%   rhow    :  ocean density (scalar variable)
%   rho     :  ice density (nodal variable)
%   g       :  gravitational acceleration
% 
%%  
        
    rho=900+zeros(MUA.Nnodes,1) ; 
    rhow=1000; 
    g=9.8/1000;
    
    
end
