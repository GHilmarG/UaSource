function [UserVar,rho,rhow,g]=DefineDensities(UserVar,CtrlVar,MUA,time,s,b,h,S,B)

%%
%  Define ice and ocean densities, as well as the gravitational acceleration.
%
% [UserVar,rho,rhow,g]=DefineDensities(UserVar,CtrlVar,MUA,time,s,b,h,S,B)
%
%  rhow    :  ocean density (scalar variable)
%  rho     :  ice density (nodal variable)
%  g       :  gravitational acceleration
% 
%%  
        
    rho=918+zeros(MUA.Nnodes,1) ; 
    rhow=1028; 
    g=9.81/1000;
    
    
end
