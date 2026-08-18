





function PlotModelAndFiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel,dudField,dvdField,dhdField,dudFieldpert,dvdFieldpert,dhdotdFieldpert)


%%
%
%
%
%
dudField=dudField(:,NodeTest);
dvdField=dvdField(:,NodeTest);
if ~isempty(dhdField)
    dhdField=dhdField(:,NodeTest);
end
%%



%% uv

UaPlots(CtrlVar,MUA,F,[dudFieldpert dvdFieldpert],FigureTitle="dv/d"+Field+" finite diff")
hold on ; plot(F.x(NodeTest)/CtrlVar.PlotXYscale,F.y(NodeTest)/CtrlVar.PlotXYscale,"o",MarkerFaceColor="r",MarkerEdgeColor="k")
title("dv/d"+Field+"finite diff. sensitivities") ; subtitle("")

UaPlots(CtrlVar,MUA,F,[dudField dvdField],FigureTitle="dv/d"+Field+" model sensitivities (exact)")
hold on ; plot(F.x(NodeTest)/CtrlVar.PlotXYscale,F.y(NodeTest)/CtrlVar.PlotXYscale,"o",MarkerFaceColor="r",MarkerEdgeColor="k")
title("dv/d"+Field+"model sensitivities (exact)") ; subtitle("")


%% dh/dt



figChdot=FindOrCreateFigure("dhdot/d"+Field+" comparision"); clf(figChdot)
T=tiledlayout("flow");

T1=nexttile;
cbar=UaPlots(CtrlVar,MUA,F,dhdField,CreateNewFigure=false)  ; title("$d\dot{h}/d"+Field+"$ sensitvity ",Interpreter="latex") ; subtitle("")
title(cbar,"")

T2=nexttile;
cbar=UaPlots(CtrlVar,MUA,F,dhdotdFieldpert,CreateNewFigure=false) ; title("$d\dot{h}/d"+Field+"$ finite differences ",Interpreter="latex") ; subtitle("") ; title(cbar,"")

if ~isempty(dhdField)
    T3=nexttile;
    UaPlots(CtrlVar,MUA,F,dhdField-dhdotdFieldpert,CreateNewFigure=false) ; title("$d\dot{h}/d"+Field+"$ differences",Interpreter="latex") ; subtitle("")
    CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);
end
T.Padding="loose";   T.TileSpacing="tight";

%%  xy-plots




figCgradu=FindOrCreateFigure("duF/d"+Field+" gradient test") ;  clf(figCgradu)
plot(dudField,dudFieldpert,"or") ;
hold on
axis equal
AX=axis;
plot([min(dudField) max(dudField)],[min(dudField) max(dudField)],"--k") ;
axis equal tight ;
xlabel(" $duF/d"+Field+"$",Interpreter="latex")  ;
ylabel("Finite difference $duF/d"+Field+"$",Interpreter="latex")
ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
axis on ; axis equal tight ; box off
title("Comparision between adjoint and finite-differences gradient calculations")
set(gcf,'Color','white')

[fitobject,gof]=fit(dudField,dudFieldpert,'poly1');
coeff=coeffvalues(fitobject);
subtitle(sprintf("slope=%f R2=%g",coeff(1),gof.rsquare))

figCgradv=FindOrCreateFigure("dvF/d"+Field+" gradient test") ;  clf(figCgradv)
plot(dvdField,dvdFieldpert,"or") ;
hold on
axis equal
plot([min(dvdField) max(dvdField)],[min(dvdField) max(dvdField)],"--k") ;
axis equal tight ;
xlabel(" $dvF/d"+Field+"$",Interpreter="latex")  ;
ylabel("Finite difference $dvF/d"+Field+"$",Interpreter="latex")
ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
axis on ; axis equal tight ; box off
title("Comparision between adjoint and finite-differences gradient calculations")
set(gcf,'Color','white')
[fitobject,gof]=fit(dvdField,dvdFieldpert,'poly1');
coeff=coeffvalues(fitobject);
subtitle(sprintf("slope=%f R2=%g",coeff(1),gof.rsquare))


if ~isempty(dhdField)

    figCgradh=FindOrCreateFigure("dhdotF/d"+Field+" gradient test") ;  clf(figCgradh)
    plot(dhdField,dhdotdFieldpert,"or") ;
    hold on
    axis equal
    plot([min(dhdField) max(dhdField)],[min(dhdField) max(dhdField)],"--k") ;
    axis equal tight ;
    xlabel(" $d\dot{h}F/d"+Field+"$",Interpreter="latex")  ;
    ylabel("Finite difference $d\dot{h}F/d"+Field+"$",Interpreter="latex")
    ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    axis on ; axis equal tight ; box off
    title("Comparision between adjoint and finite-differences gradient calculations")
    set(gcf,'Color','white')
    [fitobject,gof]=fit(dhdField,dhdotdFieldpert,'poly1');
    coeff=coeffvalues(fitobject);
    subtitle(sprintf("slope=%f R2=%g",coeff(1),gof.rsquare))
end




drawnow
input("PlotModelAndFiniteDifferenceSensitivities: Inspect, and then press RET to continue")





end