function [ab,qx,qy,dqxdx,dqxdy,dqydx,dqydy]=CalcIceShelfMeltRates(CtrlVar,MUA,u,v,s,b,S,B,rho,rhow,dsdt,as,dhdt)

%%
%
%   ab=CalcIceShelfMeltRates(CtrlVar,MUA,u,v,s,b,S,B,rho,rhow,dsdt,as,dhdt)
%
% Calculates basal meltrate from ice-flux divergence, surface mass balance (as) and surface elevation rate changes (dsdt) or alternativily from
% rates of thickness change (dhdt).
%
% Just a rough way of doing this...please look through the code and adjust as you feel needed.
%
% Calculates basal ice shelf melt rates (ab) as:
%
%  ab=dhdt+(dqxdx+dqydy)./rho-as;
%
% where:
%           dhdt is the rate of thickness change
%           dqxdx and dqydy are the flux gradients calcuated from u, v, h, and rho
%           as is the surface accumulation in meters of water equivalent
% 
% Note: It is very likely that the calculated ab values need to me smoothed or regularized. Smoothing can, for example, be done
% using the Helmholtz equation solver. To regularize one could, for example, minimize
%
%  ab' R ab + (ab-F)' P (ab-F)
%
% giving the system   (R+P) ab = P F
%
% where F is this function, and R and P regularisation and likelyhood inverse covariances. 
% 
% 
% Example would be: 
% 
%   P=MUA.M ; a=1; gs=1e7 ;  R= ga* MUA.M  + gs *(MUA.Dxx+MUA.Dyy); abEst=(R+P)\ (P*ab); 
%
% where ab was calculated by a call to this function. 
%
% Input:
%
%   u and v       : horizontal ice-velocity components
%   s, b, S and B : upper (s) and lower ice surfaces (b), ocean elevation (S) and bedrock (B)
%   rho           : ice density
%   rhow          : ocean density
%   dsdt          : rates of elevation change
%   as            : surface mass balance
%   dhdt          : rates of ice thickness change (optional, if given dsdt is ignored)
%
% Output:
%
%   ab            : lower-surface mass balance over floating areas, i.e. ice-shelf melt rate
%   qx, qy        : x and y components of ice flux
%   dqxdx, etc.   : ice-flux gradients
%
%
% Flux is calculated in units of distance/time for example as m/yr, i.e. same as velocity
%
% Sign convention: Positive values imply freezing (mass added on). Negative values imply melting (mass lost).
%
%%


if nargin==13
    fprintf('dhdt given as input. Will ignore dsdt input field\n')
else
    fprintf('calculating dhdt based on measured dsdt\n')
    dhdt=dsdt./(1-rho/rhow);
end

h=s-b;

% Flux includes density
qx=rho.*u.*h;  qy=rho.*v.*h;


% calculate flux gradients at integratino points
[dqxdx,dqxdy]=calcFEderivativesMUA(qx,MUA,CtrlVar);
[dqydx,dqydy]=calcFEderivativesMUA(qy,MUA,CtrlVar);


% Project onto nodes
dqxdx(isnan(dqxdx))=0; dqydy(isnan(dqydy))=0;
[dqxdx,dqydy]=ProjectFintOntoNodes(MUA,dqxdx,dqydy);


% And now calculate basal melt: 
ab=dhdt+(dqxdx+dqydy)./rho-as;

% Finally, put melt rates over grounded areas to zero.
GF=GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);
ab=ab.*(1-GF.node);


end