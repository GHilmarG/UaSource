function [s,b,h,S,B,rho,AGlen,n,C,m,GF,varargout]=...
    MapQuantitiesToNewFEmesh(CtrlVar,MUAnew,MUAold,hOld,time,OutsideValues,...
    varargin)



nVarargsIn = length(varargin);
fprintf('MapQuantitiesToNewFEmesh: Total number of varargin inputs = %d\n',nVarargsIn);

Nin=6 ; Nout=11;

if (nargin-Nin) ~= (nargout-Nout)
    error(' incorrect combination of input and output variables ')
end

x=MUAnew.coordinates(:,1); y=MUAnew.coordinates(:,2);


if CtrlVar.doDiagnostic
    % if a dignostic step then surface (s) and bed (b), and hence the thickness (h), are defined by the user
    [s,b,S,B,alpha]=DefineGeometry(CtrlVar.Experiment,CtrlVar,MUAnew,time,'sbSB');
    h=s-b;
else
    % if a prognostic step then surface (s) and bed (b) are defined by mapping old thickness onto
    [~,~,S,B,alpha]=DefineGeometry(CtrlVar.Experiment,CtrlVar,MUAnew,time,'SB');
    OutsideValue=0;
    h=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,x,y,OutsideValue,hOld);
end

[rho,rhow,g]=DefineDensities(CtrlVar.Experiment,CtrlVar,MUAnew,time,[],[],h,S,B);
[b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUAnew.coordinates);
GF = GL2d(B,S,h,rhow,rho,MUAnew.connectivity,CtrlVar);
[C,m]=DefineSlipperyDistribution(CtrlVar.Experiment,CtrlVar,MUAnew,time,s,b,h,S,B,rho,rhow,GF);
[AGlen,n]=DefineAGlenDistribution(CtrlVar.Experiment,CtrlVar,MUAnew,time,s,b,h,S,B,rho,rhow,GF);
%[as,ab]=DefineMassBalance(CtrlVar.Experiment,CtrlVar,MUAnew,time,s,b,h,S,B,rho,rhow,GF);

% [DTxy2,TRIxy2,DTint2,TRIint2,Xint2,Yint2,xint2,yint2,Iint2]=TriangulationNodesIntegrationPoints(MUA);

varargout=cell(nVarargsIn,1);
[varargout{:}]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,x,y,OutsideValues,varargin{:});



end