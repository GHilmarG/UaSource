function  [UserVar,s,b,S,B,rho,rhow,g]=DefineGeometryAndDensities(UserVar,CtrlVar,MUA,F,FieldsToBeDefined)


%%
%
% Defines model geometry and ice densities
%
%  [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)
%
% FieldsToBeDefined is a string indicating which return values are required. For
% example if
%
%   FieldsToBeDefined="-s-b-S-B-rho-rhow-g-"
%
% then s, b, S, B, rho, rhow and g needed to be defined.
%
% Typically, in a transient run
%
%   FieldsToBeDefined="-S-B-rho-rhow-g-"
%
% implying that only s and b do not needed to be defined, and s and b can be set to any
% value, for example s=NaN and b=NaN.
%
% It is OK to define values that are not needed, these will simply be ignored by Úa.
%
% As in all other calls:
%
%  s           is upper ice surface
%  b           is lower ice surface
%  B           is bedrock
%  S           is ocean surface
%
%   rhow    :  ocean density (scalar variable)
%   rho     :  ice density (nodal variable)
%   g       :  gravitational acceleration
%
% These fields need to be returned at the nodal coordinates. The nodal
% x and y coordinates are stored in MUA.coordinates, and also in F as F.x and F.y
%
%%


B=MismBed(F.x,F.y);

S=B*0;
b=B;
h0=300;
s=b+h0;

rho=917+zeros(MUA.Nnodes,1) ;
rhow=1030;
g=9.81/1000;



end




