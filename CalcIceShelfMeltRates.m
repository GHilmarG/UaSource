function ab=CalcIceShelfMeltRates(CtrlVar,MUA,u,v,s,b,S,B,rho,rhow,dsdt,as,dhdt)

%% Calculates basal meltrate from velocity divergence and surface mass balance and surface elevation rate changes
%
% Sign convention:
%        Positive values freezing (mass added on)
%        Negative values melting (mass lost)

if nargin==14;
    fprintf('dhdt given as input. Will ignore dsdt input field\n')
else
    dhdt=dsdt./(1-rho/rhow);
end

h=s-b;
qx=rho.*u.*h;  qy=rho.*v.*h;
[dqxdx,dqxdy]=calcFEderivativesMUA(qx,MUA,CtrlVar);
[dqydx,dqydy]=calcFEderivativesMUA(qy,MUA,CtrlVar);

dqxdx(isnan(dqxdx))=0; dqydy(isnan(dqydy))=0;
[dqxdx,dqydy]=ProjectFintOntoNodes(MUA,dqxdx,dqydy);

nOut=7;
upper=mean(dqxdx)+nOut*std(dqxdx) ; lower=mean(dqxdx)-nOut*std(dqxdx) ; 
I=dqxdx< lower | dqxdx>upper ; dqxdx(I)=NaN;

upper=mean(dqydy)+nOut*std(dqydy) ; lower=mean(dqydy)-nOut*std(dqydy) ; 
I=dqydy< lower | dqydy>upper ; dqydy(I)=NaN;



% The gradienst in fluxes tend to be very noisy, so some averaging is required.
[M,ElePerNode] = Ele2Nodes(MUA.connectivity,MUA.Nnodes);
dqxdx=Nodes2EleMean(MUA.connectivity,dqxdx);  dqxdx=M*dqxdx;
dqydy=Nodes2EleMean(MUA.connectivity,dqydy);  dqydy=M*dqydy;



%% Element averaging approach
%dqxdx=mean(dqxdx,2);  % mean values for each element
%dqydy=mean(dqydy,2);
% then project from ele to nodes by averaging over ele belonging to node
%[M,ElePerNode] = Ele2Nodes(connectivity,length(coordinates));
%dqxdxNode=M*dqxdx;
%dqydyNode=M*dqydy;
% and now smoth this over a few elements
%dqxdx=Nodes2EleMean(connectivity,dqxdxNode);  dqxdxNode=M*dqxdx;
%dqydy=Nodes2EleMean(connectivity,dqydyNode);  dqydyNode=M*dqydy;
%%

ab=dhdt+(dqxdx+dqydy)./rho-as;

% now put melt rates over grounded areas to zero.
GF=GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);
ab=ab.*(1-GF.node);


end