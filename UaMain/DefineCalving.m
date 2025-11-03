function [UserVar,LSF,c]=DefineCalving(UserVar,CtrlVar,MUA,LSF,c,F,BCs)


warning("DefineCalving:UsSource","This is the UaSource version of DefineCalving.m")


%%
%
%   [UserVar,LSF,c]=DefineCalving(UserVar,CtrlVar,MUA,LSF,c,F,BCs)
%
% Define calving the Level-Set Field (LSF) and the Calving Rate Field (c)
%
% Both the Level-Set Field (LSF) and the Calving-Rate Field (c) must be defined over the whole computational domain.
%
% The level-set option must be activated by setting
%
%  CtrlVar.LevelSetMethod=1; 
%
% in DefineInitialInputs.m
%
% The LSF should, in general, only be defined in the beginning of the run and set the initial value for the LSF. However, if
% required, the user can change LSF at any time step. The LSF is evolved by solving the Level-Set equation, so any changes
% done to LSF in this m-file will overwrite/replace the previously calculated values for LSF.
%
% The calving-rate field, c, is an input field to the Level-Set equation and needs to be defined in this m-file in each call.
%
% The variable F has F.LSF and F.c as subfields. In a transient run, these will be the corresponding values from the previous
% time step.
%
%
% In contrast to LSF, c is never evolved by Úa.  (Think of c as an input variable similar to the input as and ab for upper
% and lower surface balance, etc.)
%
% If c is returned as a NaN, ie
%
%       c=NaN;
%
% then the level-set is NOT evolved in time using by solving the level-set equation. This can be useful if, for example, the
% user simply wants to manually prescribe the calving front position at each time step.
%
%
% See more information in Ua2D_DefaultParameters.
%
%%

%% initialize LSF
if isempty(F.LSF)   % Do I need to initialize the level set function?

    % If, for example, all floating ice should be removed at the beginning of the simulations
    % then the level-set can simply be defined in such a way that LSF>0 where the ice is grounded, and <0 otherwise

    if contains(UserVar.RunType,"-c0isGL0-")  % -> Initial calving front (c0) is set a initial grounding line position (GL0)

        LSF=-ones(MUA.Nnodes,1) ;
        LSF(F.GF.node>0.5)=+1;    % LSF set dependent on the values grounded-floating mask
        Xc=[] ;  % If Xc and Yc are left empty, the Xc and Yc will be calculated as the zero control of the LSF field
        Yc=[] ; 

    else

        Xc=UserVar.CalvingFront0.Xc;  % Here is is assumed that the user has defined the desired initial location of the calving front, and stored those in UserVar
        Yc=UserVar.CalvingFront0.Yc;

        % A rough sign-correct initialization for the LSF
        io=inpoly2([F.x F.y],[Xc(:) Yc(:)]);
        LSF=-ones(MUA.Nnodes,1) ;
        LSF(io)=+1;

    end

    % figure ; PlotMuaMesh(CtrlVar,MUA);   hold on ; plot(F.x(io)/1000,F.y(io)/1000,'or')

    % Not really needed, but it might then be good to initialize the level set so that LSF is approximately a signed distance
    % function.  
    %
    [xc,yc,LSF]=CalvingFrontLevelSetGeometricalInitialisation(CtrlVar,MUA,Xc,Yc,LSF,plot=true,ResampleCalvingFront=true);


end

%% Define calving rate (if needed)

if  CtrlVar.LevelSetEvolution=="-Prescribed-"

    c=nan;   % setting the calving rate to nan implies that the level set is not evolved

elseif  CtrlVar.CalvingLaw.Evaluation=="-int-"

    % Generally, one would define the calving rate at the nodes, unless the calving rate is a function of the level-set itself.
    %
    % It the calving rate, c, depends on the level set function, LSF, then various derivatives involving c and LSF (phi) need to be
    % calculated as well for the Newton-Raphson methods to achieve second-order convergence

    % If the calving law is defined at nodes, you still need to use DefineCalving.m as well to define LSF. 

    % Generally there is no need to call 'DefineCalvingAtIntegrationPoints' from within DefineCalving.m and the calving rate, c,
    % returned by DefineCalving.m is not used.
    %
    % However, one might find it convenient to define c here as otherwise c is not defined at nodes and can not be plotted at a
    % nodal variable.
    %
    % 
    %
    c=DefineCalvingAtIntegrationPoints(UserVar,CtrlVar,nan,nan,F) ;

else

    % Here the calving rate is defined at the nodes. This is presumably the most typical case
   
    CliffHeight=min((F.s-F.S),F.h) ;
    c=10*CliffHeight ;  % an example of calving law which depends linearly on the freeboard height. 


end

