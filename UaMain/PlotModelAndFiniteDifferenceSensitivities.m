





function PlotModelAndFiniteDifferenceSensitivities(CtrlVar,MUA,BCs,F,l,Field,NodeTest,dudField,dvdField,dhdField,dudFieldpert,dvdFieldpert,dhdotdFieldpert,SubtitleString)

narginchk(14,14)

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
title(sprintf("(du/d"+Field+",dv/d"+Field+") finite diff. sensitivities at node %i",NodeTest)) ; 
subtitle(SubtitleString,Interpreter="latex")

UaPlots(CtrlVar,MUA,F,[dudField dvdField],FigureTitle="dv/d"+Field+" model sensitivities (exact)")
hold on ; plot(F.x(NodeTest)/CtrlVar.PlotXYscale,F.y(NodeTest)/CtrlVar.PlotXYscale,"o",MarkerFaceColor="r",MarkerEdgeColor="k")
title(sprintf("(du/d"+Field+",dv/d"+Field+") model sensitivities at node %i",NodeTest)) ; 
subtitle(SubtitleString,Interpreter="latex")


%% dh/dt


figChdot=FindOrCreateFigure("Test: dhdot/d"+Field+" comparision"); clf(figChdot)
T=tiledlayout("flow");

T1=nexttile;
cbar=UaPlots(CtrlVar,MUA,F,dhdField,CreateNewFigure=false)  ; 
title("$d\dot{h}/d"+Field+"$ sensitvity ",Interpreter="latex") ; 
subtitle(SubtitleString,interpreter="latex")
title(cbar,"")

T2=nexttile;
cbar=UaPlots(CtrlVar,MUA,F,dhdotdFieldpert,CreateNewFigure=false) ; 
title("$d\dot{h}/d"+Field+"$ finite differences ",Interpreter="latex") ;  title(cbar,"")
subtitle(SubtitleString,interpreter="latex")

if ~isempty(dhdField)
    T3=nexttile;
    UaPlots(CtrlVar,MUA,F,dhdField-dhdotdFieldpert,CreateNewFigure=false) ; 
    title("Test: $d\dot{h}/d"+Field+"$ differences",Interpreter="latex") ; 
    subtitle(SubtitleString,interpreter="latex")
    CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);
end
T.Padding="loose";   T.TileSpacing="tight";

%%  xy-plots




figCgradu=FindOrCreateFigure("Test: duF/d"+Field+" sensitivity") ;  clf(figCgradu)
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
title("Comparision between adjoint and finite-differences sensitivities")
set(gcf,'Color','white')

[fitobject,gof]=fit(dudField,dudFieldpert,'poly1');
coeff=coeffvalues(fitobject);
subtitle(SubtitleString+sprintf("  slope=%f $R^2$=%g",coeff(1),gof.rsquare),Interpreter="latex")

Diff_u=norm(dudField-dudFieldpert)/norm(dudField);
fprintf("%s: normalized differences in u sensitvities is: %g\n",Field,Diff_u)    






figCgradv=FindOrCreateFigure("Test: dvF/d"+Field+" sensitiviy") ;  clf(figCgradv)
plot(dvdField,dvdFieldpert,"or") ;
hold on
axis equal
plot([min(dvdField) max(dvdField)],[min(dvdField) max(dvdField)],"--k") ;
axis equal tight ;
xlabel(" $dvF/d"+Field+"$",Interpreter="latex")  ;
ylabel("Finite difference $dvF/d"+Field+"$",Interpreter="latex")
ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
axis on ; axis equal tight ; box off
title("Comparision between adjoint and finite-differences sensitvities")
set(gcf,'Color','white')
[fitobject,gof]=fit(dvdField,dvdFieldpert,'poly1');
coeff=coeffvalues(fitobject);
subtitle(SubtitleString+sprintf("  slope=%f $R^2$=%g",coeff(1),gof.rsquare),Interpreter="latex")

Diff_v=norm(dvdField-dvdFieldpert)/norm(dvdField);
fprintf("%s: normalized differences in v sensitvities is: %g\n",Field,Diff_v)    


if ~isempty(dhdField)

    figCgradh=FindOrCreateFigure("Test: dhdotF/d"+Field+" sensitivity") ;  clf(figCgradh)
    plot(dhdField,dhdotdFieldpert,"or") ;
    hold on
    axis equal
    plot([min(dhdField) max(dhdField)],[min(dhdField) max(dhdField)],"--k") ;
    axis equal tight ;
    xlabel(" $d\dot{h}F/d"+Field+"$",Interpreter="latex")  ;
    ylabel("Finite difference $d\dot{h}F/d"+Field+"$",Interpreter="latex")
    ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    axis on ; axis equal tight ; box off
    title("Comparision between adjoint and finite-differences sensitivities")
    set(gcf,'Color','white')
    [fitobject,gof]=fit(dhdField,dhdotdFieldpert,'poly1');
    coeff=coeffvalues(fitobject);
    subtitle(SubtitleString+sprintf("  slope=%f $R^2$=%g",coeff(1),gof.rsquare),Interpreter="latex")


    Diff_hdot=norm(dhdField-dhdotdFieldpert)/norm(dhdField);
    fprintf("%s: normalized differences in dh/dt sensitvities is: %g\n",Field,Diff_hdot)    

end




drawnow

% fprintf("PlotModelAndFiniteDifferenceSensitivities: Inspect in debugger and then continue: [F5] \n")
% keyboard
% 






end