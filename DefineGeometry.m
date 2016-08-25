
function [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)

%%
%  Defines geometrical variables
% [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)
%
% Outputs:
%
% * s : upper glacier surface 
% * b : lower glacier surface 
% * S : ocean surface 
% * B : bedrock
% * alpha : tilt of the coordinate system with respect to gravity.
% 
% s, b, S and B are nodal variables.
%
%
% FieldsToBeDefined    : a string indicating which output variables are required
%                        at the current stage of the run.
%
% If for example this string is 'sbSB' , then s, b, S and B are required. If the
% string is 'B'  then only B is required. The other variables must be : defined
% as well but their values will not be used.
%
%                      
% Note: Floating conditions are taken care of internally and s and b will be
% modified accordingly while conserving thickness.
%
% Note: If for some reason s and b are defined so that s-b<CtrlVar.ThickMin then
% s and b will be modified internally to ensure that min ice thickness is no
% less than CtrlVar.ThickMin
%
% Note: In a transient run s and b are updated internally and do not need to be
% defined by the user (in this case FieldsToBeDefined='SB') except at the start
% of the run.
%
% Note: After remeshing, variables are generally mapped internally from the old
% to the new mesh using some interpolation methods (Usually a natural neighbour
% scattered interpolant). But in a diagnostic/static run, geometric variables
% (s,b,S,B) are always defined over the new mesh through a call to
% [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined).
%       
%%
    x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
    
    alpha=0.01; % tilt of the coordinate system with respect to gravity
    
    
    hmean=1000; 
    ampl_b=0.5*hmean; sigma_bx=5000 ; sigma_by=5000;
    Deltab=ampl_b*exp(-((x/sigma_bx).^2+(y/sigma_by).^2));
    Deltab=Deltab-mean(Deltab);
    
    B=zeros(MUA.Nnodes,1) + Deltab;  % Gaussian bedrock
    S=B*0-1e10;                      % Put ocean surface so low that no sections of the glacier will be afloat
    b=B;                             % Set lower glacier surface elevation to that of bedrock
    s=B*0+hmean;                     % Upper glacier surface.
    
    
end
