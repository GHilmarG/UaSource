function [UserVar,uo,vo,Co,mo,ua,va,Ca,ma]=DefineSeaIceParameters(UserVar,CtrlVar,MUA,BCs,GF,ub,vb,ud,vd,uo,vo,ua,va,s,b,h,S,B,rho,rhow,AGlen,n,C,m)

%%
%
% This m-file is used to define drag due to prescribed ocean current and wind speed.
%
%
% 
% The relationship between (precribed) currents and winds and the resulting drag is formally identical to Weertman sliding law, but
% only acts over the floating sections of the domain.
%
%
%   (uo,vo)              ocean velocity, defined on nodes, can be spatially variable
%   (ua,va)              wind velocity, defined on nodes, can be spatially variable
%
%   Co , mo              parameters related to ocean-induced drag.  See the UaCompendium for exact defintion
%   Ca , ma              parameters related to wind-induced drag.  See the UaCompendium for exact defintion
%
%
% For example:
% 
%   uo=zeros(MUA.Nnodes,1);
%   vo=zeros(MUA.Nnodes,1);
%   ua=zeros(MUA.Nnodes,1);
%   va=zeros(MUA.Nnodes,1);
% 
%   mo=2; Co=1e+10; 
%   ma=1; Ca=1e+10;
%
%
%
%%


Co=[];
mo=[];

Ca=[];
ma=[];




end