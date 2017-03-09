function [UserVar,EleSizeDesired,ElementsToBeRefined]=...
    DefineDesiredEleSize(UserVar,CtrlVar,MUA,x,y,EleSizeDesired,ElementsToBeRefined,s,b,S,B,rho,rhow,ub,vb,ud,vd,GF,NodalErrorIndicators)

%%
% Define desired sizes of elements or specify which elements to refine.
%
% [EleSizeDesired,ElementsToBeRefined]=DefineDesiredEleSize(CtrlVar,MUA,x,y,EleSizeDesired,ElementsToBeRefined,s,b,S,B,rho,rhow,ub,vb,ud,vd,GF,NodalErrorIndicators)
%
% Only used in combination with adaptive meshing.
%
% Allows user to set desired ele sizes (EleSize) at given locations (x,y).
%
% On input x, y, EleSize are desired ele sizes at (x,y) as
% calculated by Úa based on some user-defined criteria.
%
% On output x,y,EleSize are user-modified values.
%
% Do not modify the size of the (nodal) vector `EleSizeDesired' or the logical (element)
% vector 'ElementsToBeRefine', only the values.
%
% x,y are the locations where new element sizes are specifed. 
%
% ElementsToBeRefined can either be a logical array in which case values set to true/1 indicate elements
% to be refined, or a list of numbers of elements to be refined.
%%
 
 
    
end
