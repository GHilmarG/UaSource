function [ab,qx,qy,dqxdx,dqxdy,dqydx,dqydy]=CalcIceShelfMeltRates(CtrlVar,MUA,u,v,s,b,S,B,rho,rhow,dsdt,as,dhdt)

%%
% ab=CalcIceShelfMeltRates(CtrlVar,MUA,u,v,s,b,S,B,rho,rhow,dsdt,as,dhdt)
% Calculates basal meltrate from velocity divergence and surface mass balance and surface elevation rate changes
%
% Just a rough way of doing this...
%
%  Calculates basal ice shelf melt rates (ab) as:
%
%  ab=dhdt+(dqxdx+dqydy)./rho-as;
%
%  where:
%           dhdt is the rate of thickness change
%           dqxdx and dqydy are the flux gradients calcuated from u, v, h, and rho
%           as is the surface accumulation in meters of water equivalent
%
%           Flux is calculated in units of distance/time for example as m/yr, i.e. same as velocity
%
% Sign convention:
%        Positive values freezing (mass added on)
%        Negative values melting (mass lost)

if nargin==13;
    fprintf('dhdt given as input. Will ignore dsdt input field\n')
else
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



% some rather ad-hoc eliminayion of outliers
% through out data more than nErr standard deviation from mean
nErr=7;
upper=mean(dqxdx)+nErr*std(dqxdx) ; lower=mean(dqxdx)-nErr*std(dqxdx) ; 
I=dqxdx< lower | dqxdx>upper ; dqxdx(I)=NaN;

upper=mean(dqydy)+nErr*std(dqydy) ; lower=mean(dqydy)-nErr*std(dqydy) ; 
I=dqydy< lower | dqydy>upper ; dqydy(I)=NaN;



% The gradienst in fluxes tend to be very noisy, so some averaging is required.
% Again this is rather ad-hoc
[M,ElePerNode] = Ele2Nodes(MUA.connectivity,MUA.Nnodes);
dqxdx=Nodes2EleMean(MUA.connectivity,dqxdx);  dqxdx=M*dqxdx;
dqydy=Nodes2EleMean(MUA.connectivity,dqydy);  dqydy=M*dqydy;



ab=dhdt+(dqxdx+dqydy)./rho-as;

% now put melt rates over grounded areas to zero.
GF=GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);
ab=ab.*(1-GF.node);


end