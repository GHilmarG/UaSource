





function [UserVar,RunInfo,F,l,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=...
    DefineDesiredEleSize(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,x,y,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,NodalErrorIndicators)


%%
% Define desired sizes of elements or specify which elements to refine or coarsen.
%
%   [UserVar,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=...
%            DefineDesiredEleSize(UserVar,CtrlVar,MUA,x,y,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,s,b,S,B,rho,rhow,ub,vb,ud,vd,GF,NodalErrorIndicators)
%
% Only used in combination with adaptive meshing.
%
% You need to set
%
%   CtrlVar.AdaptMesh=1;  
%
% in Ua2D_InitialUserinput for this m-file to be called.
%
% Allows user to set:
% 
% * EleSizeDesired when using global mesh refinement
% * ElementsToBeRefined when using local mesh refinement with either the red-green or the newest vertex bisection
% * ElementsToBeRefined and ElementsToBeCoarsened when using local mesh refinement with the the newest vertex bisection
% 
%
% On input EleSize are desired element sizes at (x,y) as
% calculated by Úa based on some user-defined criteria.
%
% On output EleSize are user-modified values.
%
% Do not modify the size of the (nodal) vector `EleSizeDesired' or the logical (element)
% vector 'ElementsToBeRefine', only the values.
%
% When using the global remeshing option x,y are the locations where new element sizes are specified (these are the coordinates of the mesh)
% 
% *Note: When using the local remeshing option, x and y as given on input are not relevant. 
%       In this case use MUA.xEle and MUA.yEle as the x, y locations where the elements are to be refined or coarsened.* 
%
% ElementsToBeRefined can either be a logical array in which case values set to true/1 indicate elements
% to be refined, or a list of numbers of elements to be refined.
%
% Note that this m-file is only called if the adaptive meshing option is used.
% Also, that elements will only be refined/coarsened if local mesh refinement is
% used. These options must be set accordingly in Ua2D_InitialUserInput.
%
% 
% *Example:* To set desired element size to 1000 within a given boundary (this boundary
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
% Here Boundary does not have to be just a simple square, it can be a polygon of any shape.   
%
% *Example:* To set all ele size of all floating elements (i.e. ice shelves)
% to 1000:
%
%   EleSizeDesired(GF.Node<0.5)=1000;
%
% *Example* defining either EleSizeDesired or ElementsToBeRefined depending on adaptive meshing option selected:
%
%     switch lower(CtrlVar.MeshRefinementMethod)
%     
%         case 'explicit:global' 
%         
%         % when using global mesh refinement, return EleSizeIndicator defined at nodes
% 
%             EleSizeIndicator=EleSizeDesired;
% 
%             EleSizeIndicator(GF.node<0.1)=UserVar.MeshSizeIceShelves;
%             EleSizeDesired=min(EleSizeDesired,EleSizeIndicator);
%         
%             EleSizeIndicator(s<1500)=CtrlVar.MeshSizeMax/5;
%             EleSizeDesired=min(EleSizeDesired,EleSizeIndicator);
%         
%             xmin=-1727e3   ; xmax=-1100e3 ; ymin=-600e3 ; ymax=-20.e3;
%             ind=x< xmax & x>xmin & y>ymin & y< ymax ;
%             EleSizeDesired(~ind)=CtrlVar.MeshSizeMax;
%         
%         case 'explicit:local:newest vertex bisection'
%         
%         % When using local mesh refinement, return ElementsToBeRefined and ElementsToBeCoarsened defined over elements
%         %
%         % ElementsToBeCoarsened is only used in combination with the 'newest vertex bisection' local mesh-refinement method 
%         %
%             xmin=-1727e3   ; xmax=-1100e3 ; ymin=-600e3 ; ymax=-20.e3;
%             ind=MUA.xEle < xmax & MUA.xEle > xmin & MUA.yEle >ymin & MUA.yEle < ymax ;
%       
%             ElementsToBeRefined(~ind)=false; 
% 
%     end
% 
% 
% 
%%
 
 
    
end
