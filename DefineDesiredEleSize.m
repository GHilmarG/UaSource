function [UserVar,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=...
            DefineDesiredEleSize(UserVar,CtrlVar,MUA,x,y,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,s,b,S,B,rho,rhow,ub,vb,ud,vd,GF,NodalErrorIndicators)
        

%%
% Define desired sizes of elements or specify which elements to refine or
% coarsen.
%
%   [EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=...
%       DefineDesiredEleSize(CtrlVar,MUA,x,y,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,s,b,S,B,rho,rhow,ub,vb,ud,vd,GF,NodalErrorIndicators)
%
% Only used in combination with adaptive meshing.
%
% Allows user to set desired ele sizes (EleSize) at given locations (x,y).
%
% On input EleSize are desired ele sizes at (x,y) as
% calculated by Úa based on some user-defined criteria.
%
% On output EleSize are user-modified values.
%
% Do not modify the size of the (nodal) vector `EleSizeDesired' or the logical (element)
% vector 'ElementsToBeRefine', only the values.
%
% x,y are the locations where new element sizes are specifed. 
%
% ElementsToBeRefined can either be a logical array in which case values set to true/1 indicate elements
% to be refined, or a list of numbers of elements to be refined.
%
% Note that this m-file is only called if the adaptive meshing option is used.
% Also, that elements will only be refined/coarsened if local mesh refinement is
% used. These options must be set accordingly in Ua2D_InitialUserInput.
%
% 
% *Example:* To set desired ele sized to 1000 within a given boundary (this boundary
% must of course be within the overall boundary of the computational
% domain):
%
%   Boundary=[0        0 ; ...
%           10e3      0 ; ...
%           10e3      10e3;
%           0       10e3];
% 
%   I=inpoly([x y],Boundary) ;
%   EleSizeDesired(I)=1000; 
%
% Here Boundary doese not have to be just a simple square, it can be a polygon of any shape.   
%
%
% 
%
% *Example:* To set all ele size of all floating elements (i.e. ice shelves)
% to 1000:
%
%   EleSizeDesired(GF.Node<0.5)=1000;
%
%%
 
 
    
end
