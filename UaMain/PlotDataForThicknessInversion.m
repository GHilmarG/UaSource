











function PlotDataForThicknessInversion(UserVar,CtrlVar,MUA,F,BCs,Priors,Meas,htrue) 


narginchk(8,8)



cbar=UaPlots(CtrlVar,MUA,F,Meas.h,FigureTitle="h measured") ;
xlabel("$x$ (km)  ") ; ylabel("$y$ (km)  ") ; title(cbar,"(m a.s.l.)")
title("Measured Ice thickness",Interpreter="latex")
hold on  ; PlotGroundingLines(); 

cbar=UaPlots(CtrlVar,MUA,F,"-uv-",FigureTitle="measured velocities");
xlabel("$x$ (km)  ") ; ylabel("$y$ (km)  ")
title("Measured surface velocities",Interpreter="latex")
hold on  ; PlotGroundingLines(); 

cbar=UaPlots(CtrlVar,MUA,F,F.as,FigureTitle="measured mass balance");
xlabel("$x$ (km)  ") ; ylabel("$y$ (km)  ") ; title(cbar,"(m/yr)",interpreter="latex")
title("Surface mass balance",Interpreter="latex")
hold on  ; PlotGroundingLines(); 

cbar=UaPlots(CtrlVar,MUA,F,F.dsdt,FigureTitle="surface rate of elevation changes");
title("Surface elevation changes",Interpreter="latex")
title(cbar,"(m/yr)",interpreter="latex")
xlabel("$x$ (km)  ") ; ylabel("$y$ (km)  ")
hold on  ; PlotGroundingLines(); 













end