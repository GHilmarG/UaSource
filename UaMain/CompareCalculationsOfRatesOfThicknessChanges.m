





function CompareCalculationsOfRatesOfThicknessChanges(UserVar,RunInfo,CtrlVar,MUA,F,F0,BCs,l)

%% Calculates thickness changes by:
%
%
% 
% 
% # estimating dh/dt directly from evaluation the mass-conservation based on current geometry alone.
% # estimating dh/dt by solving the mass-conservation equation to time-step from t=t0 to t=t, using the theta method.
% 
%
%
%%


%% Conclusions:
%
% Based on a test done with the Thule geometry, all three produce almost identical results when using SUPG.tau="taut","tau1"
% and "tau2"; Explicit agrees less well, but this is almost certainly because the explicit method does not use the level set
% approach, nor the implicit mass-balance feedback related to thickness min. and level-set.
%
% There are some,small, but numerically significant differences in dh/dt estimates , even between the -h- solve and the -uvh-
% solve. These are primarily around the grounding line, and I surmise that this is do to inherent differences in the
% treatment of the grounding-line motion (no movement in the -h- solve, movement in the -uvh- solve) ), but this needs to be
% tested further.
%
%
%%


%% Explicit estimate based on solving the mass-conservation equation for dh/dt
[UserVar,dhdtExplicit]=dhdtExplicitSUPG(UserVar,CtrlVar,MUA,F,BCs);

%% Implicit estimate based on solving the mass-conservation equation implicitly for h
CtrlVar.ThicknessConstraints=true;
CtrlVar.LevelSetMethodAutomaticallyApplyMassBalanceFeedback=true;
CtrlVar.ThicknessPenalty=true;
% CtrlVar.h.SUPG.tau="tau0";
% CtrlVar.h.SUPG.tau="tau0";
% CtrlVar.h.SUPG.tau="tau2";  % smooth 
% CtrlVar.h.SUPG.tau="tau1";  % smooth
% CtrlVar.h.SUPG.tau="taus";  % wiggles 
% CtrlVar.h.SUPG.tau="taut";  % smooth

[UserVar,RunInfo,h1]=MassContinuityEquationNewtonRaphsonThicknessContraints(UserVar,RunInfo,CtrlVar,MUA,F0,F,l,BCs) ;

dhdtImplicit=(h1-F0.h)/F.dt;

%%

F.dhdt=(F.h-F0.h)/F.dt;

% %% Limit to grounded ice?
% dhdtExplicit=dhdtExplicit.*F.GF.node;
% dhdtImplicit=dhdtImplicit.*F.GF.node;
% F.dhdt=F.dhdt.*F.GF.node;
% 
% %% Limit to interior?
% R=sqrt(F.x.*F.x+F.y.*F.y) ;
% I=R<50e3;
% dhdtExplicit(~I)=0;
% dhdtImplicit(~I)=0;
% F.dhdt(~I)=0;
% %% Plots


UaPlots(CtrlVar,MUA,F,dhdtExplicit,FigureTitle="dhdt explicit",GetRidOfValuesDownStreamOfCalvingFronts=false)
CL=clim; 
CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);

UaPlots(CtrlVar,MUA,F0,dhdtImplicit,FigureTitle="dhdt implicit",GetRidOfValuesDownStreamOfCalvingFronts=false)
plot(F.x(BCs.hPosNode)/1000,F.y(BCs.hPosNode)/1000,"or",MarkerSize=2)
clim(CL)
CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);


UaPlots(CtrlVar,MUA,F,F.dhdt,FigureTitle="F.dhdt",GetRidOfValuesDownStreamOfCalvingFronts=false)
clim(CL)
plot(F.x(BCs.hPosNode)/1000,F.y(BCs.hPosNode)/1000,"or",MarkerSize=2)
CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
title("dh/dt as in F.dhdt")

DeltaExplicit=dhdtExplicit-F.dhdt;
UaPlots(CtrlVar,MUA,F,DeltaExplicit,FigureTitle="dhdtExplicit-F.dhdt",GetRidOfValuesDownStreamOfCalvingFronts=false)
clim(CL)
CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
title("dh/dt explicit - dh/dt as in F.dhdt")

DeltaImplicit=dhdtImplicit-F.dhdt;
UaPlots(CtrlVar,MUA,F,DeltaImplicit,FigureTitle="dhdtMass-F.dhdt",GetRidOfValuesDownStreamOfCalvingFronts=false)
clim(CL)
CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
title("dh/dt mass - dh/dt as in F.dhdt")

%% Norm of difference


DE= sqrt(DeltaExplicit'*MUA.M* DeltaExplicit)/MUA.Area; 
DI= sqrt(DeltaImplicit'*MUA.M* DeltaImplicit)/MUA.Area; 

fprintf("DE=%f \t DM=%f \n",DE,DI)

FindOrCreateFigure("DE/DM")
hold on
plot(F.time,DE,"r",Marker="square",LineStyle="none")
plot(F.time,DI,"b",Marker="diamond",LineStyle="none")
xlabel("time")
ylabel("Explicit diff (red) / Implicit diff (blue)")

end