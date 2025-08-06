


function [UserVar,ElementsToBeDeactivated]=DefineElementsToDeactivate(UserVar,RunInfo,CtrlVar,MUA,F,BCs,ElementsToBeDeactivated) 

%%  Manually deactivate elements within a mesh.   
%
% This file is called within each run step if CtrlVar.ManuallyDeactivateElements=true
%  
%   ElementsToBeDeactivated  : a list of elements to be deactivated.
%                              this can either be a logical variable or a vector with element numbers that are to be
%                              deactivated.
% 
% *Example:* To deactivate elements 10 and 23, set: 
%
%   ElementsToBeDeactivated=false(MUA.Nele,1) ; 
%   ElementsToBeDeactivated([10;23])=true ;
%
% *Example:* To deactivate elements with element x coordinates larger than 1000, set: 
%
%  ElementsToBeDeactivated=MUA.xEle > 1000
% 
%
% 
%
% *Example:*  To deactivate all elements outside of the region bounded by
% BoundaryCoordinates in run-step 2:
%
% 
%   if CtrlVar.CurrentRunStepNumber==2
%
%       BoundaryCoordinates=[0 0 ; 10e3 0 ; 10e3 20e3 ; -5e3 15e3 ] ;    
%       In=inpoly([MUA.xEle MUA.yEle],BoundaryCoordinates);
%       ElementsToBeDeactivated=~In;
%     
%     
%       figure
%       PlotMuaMesh(CtrlVar,MUA,ElementsToBeDeactivated,'r')
%       hold on
%       PlotMuaMesh(CtrlVar,MUA,~ElementsToBeDeactivated,'k')
%  
%   end
% 
% *Example:* To deactivate elements where none of the nodes have ice thickness greater
% than 2*CtrlVar.ThickMin : 
%
%   AboveMinThickNodes = F.h > 2*CtrlVar.ThickMin ; 
%   MinThickElements=AllElementsContainingGivenNodes(MUA.connectivity,AboveMinThickNodes) ; 
%   NewElementsToBeDeactivated=~MinThickElements ; 
%  ElementsToBeDeactivated=ElementsToBeDeactivated | NewElementsToBeDeactivated ; % Here it is assumed that this is a logical list
%

% *Example:* To plot elements that are to be deactivated: 
% 
%   UaPlots(CtrlVar,MUA,F,F.h,GetRidOfValuesDownStreamOfCalvingFronts=false,FigureTitle="Elements to be deactivated"); 
%   set(gca,'ColorScale','log') 
%   hold on 
%   PlotMuaMesh(CtrlVar,MUA);
%   hold on
%   PlotMuaMesh(CtrlVar,MUA,ElementsToBeDeactivated,color="r",LineWidth=2)
%   title("Elements to be deactivated in red")
% 
%
% Note that the x and y coordinates of the elements, defined as the mean x and y values of nodes, can be found in :
%
%   MUA.xEle
%   MUA.yEle
%
% 
%%




end