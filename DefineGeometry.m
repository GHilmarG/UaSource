

function [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)

%%
%   Defines model geometry
%
%  [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)
%
% FieldsToBeDefined is a string indicating which return values are required. For
% example if
%
%   FieldsToBeDefined='sbSB'
%
% then s, b, S and B needed to be defined.
%
% Typically, in a transient run
%
%   FieldsToBeDefined='SB'
%
% implying that only S and B needed to be defined, and s and b can be set to any
% value, for example s=NaN and b=NaN.
%
% As in all other calls:
%
%  s           is upper ice surface 
%  b           is lower ice surface 
%  B           is bedrock
%  S           is ocean surface
%
% These fields need to be returned at the nodal coordinates. The nodal
% coordinates are stored in MUA.coordinates
%
% alpha         is the tilt of the coordinate system with respect to gravity
%               (not the slope of the ice surface). alpha
%               is a scalar variable, and usully alpha=0
%
%%

x=MUA.coordinates(:,1); 
y=MUA.coordinates(:,2);
alpha=0.;


B=MismBed(x,y);

S=B*0;
b=B;
h0=300;
s=b+h0;

end
