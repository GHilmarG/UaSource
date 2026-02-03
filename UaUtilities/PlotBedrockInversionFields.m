



function PlotBedrockInversionFields(CtrlVar,MUA,F,Priors,InvFinalValues,InvStartValues,Meas)

narginchk(7,7)


if ~isempty(Priors.TrueB)

    % the True B field is not empty. This could, for example, be a synthetic test case.

    figB=FindOrCreateFigure("True and estimated B"); clf(figB)
    TB=tiledlayout("flow") ;

    ax1=nexttile;
    UaPlots(CtrlVar,MUA,F,Priors.TrueB,CreateNewFigure=false);
    hold on ; PlotGroundingLines(); PlotCalvingFronts();
    title("True B")
    subtitle("")

    ax2=nexttile;
    UaPlots(CtrlVar,MUA,F,InvFinalValues.B,CreateNewFigure=false);
    hold on ; PlotGroundingLines(); PlotCalvingFronts();
    title("Retrieved B")
    subtitle("")


    if isempty(Priors.TrueB)
        Priors.TrueB=F.b+nan;
    end


    ax3=nexttile;
    UaPlots(CtrlVar,MUA,F,InvFinalValues.B-Priors.TrueB,CreateNewFigure=false);
    hold on ; PlotGroundingLines(); PlotCalvingFronts();
    title("B retrieved - B true")
    subtitle("")


    ax4=nexttile;
    UaPlots(CtrlVar,MUA,F,InvStartValues.B,CreateNewFigure=false);
    hold on ; PlotGroundingLines(); PlotCalvingFronts();
    title("B at start of current inversion run")
    subtitle("")

    ax5=nexttile;
    UaPlots(CtrlVar,MUA,F,Priors.B,CreateNewFigure=false);
    hold on ; PlotGroundingLines(); PlotCalvingFronts();
    title("B prior") ; subtitle("")

    ax6=nexttile;
    UaPlots(CtrlVar,MUA,F,InvFinalValues.B-Priors.B,CreateNewFigure=false);
    hold on ; PlotGroundingLines(); PlotCalvingFronts();
    title("Retrieved B -  Prior B ")
    subtitle("")


    % feels a bit clumsy way of ensuring that each tile as its own colorbar, but I can't think of a simpler approach
    set(figB,CurrentAxes=ax1) ;
    cl=clim;
    if min(cl) <0 && max(cl)> 0
        CM=cmocean('balanced',25,'pivot',0) ; colormap(ax1,CM);
    else
        CM=cmocean('balanced',25) ; colormap(ax1,CM);
    end

    set(figB,CurrentAxes=ax2) ;
    cl=clim;
    if min(cl) <0 && max(cl)> 0
        CM=cmocean('balanced',25,'pivot',0) ; colormap(ax2,CM);
    else
        CM=cmocean('balanced',25) ; colormap(ax2,CM);
    end

    set(figB,CurrentAxes=ax3) ;

    cl=clim;
    if min(cl) <0 && max(cl)> 0
        CM=cmocean('balanced',25,'pivot',0) ; colormap(ax3,CM);
    else
        CM=cmocean('balanced') ; colormap(ax3,CM);
    end

    set(figB,CurrentAxes=ax4) ;
    cl=clim;
    if min(cl) <0 && max(cl)> 0
        CM=cmocean('balanced',25,'pivot',0) ; colormap(ax4,CM);
    else
        set(figB,CurrentAxes=ax4) ;  CM=cmocean('balanced',25) ; colormap(ax4,CM);
    end

    set(figB,CurrentAxes=ax5) ;

    cl=clim;
    if min(cl) <0 && max(cl)> 0
        CM=cmocean('balanced',25,'pivot',0) ; colormap(ax5,CM);
    else
        CM=cmocean('balanced',25) ; colormap(ax5,CM);
    end

    set(figB,CurrentAxes=ax6) ;
    cl=clim;
    if min(cl) <0 && max(cl)> 0
        CM=cmocean('balanced',25,'pivot',0) ; colormap(ax6,CM);
    else
        CM=cmocean('balanced',25) ; colormap(ax6,CM);
    end


    % figB.Position=[200 200 1300 800];
    TB.TileSpacing="tight";
    TB.Padding="tight";
    %colormap(othercolor("Mdarkterrain",32))

end


if ~isempty(Meas.B)

    % the Meas.B field is not empty. These are measurements projected onto nodes

    BDiff=InvFinalValues.B-Meas.B ;
    BErr=sqrt(spdiags(Meas.BCov));
    I=isfinite(BErr);


    figB=FindOrCreateFigure("Measured and estimated B"); clf(figB)
    TB=tiledlayout("flow") ;

    ax1=nexttile;
    UaPlots(CtrlVar,MUA,F,Meas.B,CreateNewFigure=false);
    climMeasured=clim;
    %hold on ; PlotGroundingLines(); PlotCalvingFronts();
    title("Measured B (projected onto nodes)")
    subtitle("")
    hold on
    plot(F.x(I)/CtrlVar.PlotXYscale,F.y(I)/CtrlVar.PlotXYscale,".k",MarkerSize=3,DisplayName="Meas. locations")
    lg=legend;
    lg.String{1}="$(B_{\mathrm{Retrieved}}$";
    lg.Interpreter="latex";

    ax2=nexttile;
    UaPlots(CtrlVar,MUA,F,InvFinalValues.B,CreateNewFigure=false);
    %hold on ; PlotGroundingLines(); PlotCalvingFronts();
    title("Retrieved B") ; subtitle("")
    clim(climMeasured);


    ax3=nexttile;
    UaPlots(CtrlVar,MUA,F,BDiff,CreateNewFigure=false);
    hold on
    plot(F.x(I)/CtrlVar.PlotXYscale,F.y(I)/CtrlVar.PlotXYscale,".k",MarkerSize=3,DisplayName="Meas. locations")
    %hold on ; PlotGroundingLines(); PlotCalvingFronts();
    title("B estimated - B measured") ; subtitle("")
    lg=legend;
    lg.String{1}="$(B_{\mathrm{Retrieved}}-B_{\mathrm{Meas}})$";
    lg.Interpreter="latex";

    BErr(~I)=nan;
    BDiffErr=BDiff./BErr   ;

    ax4=nexttile;
    UaPlots(CtrlVar,MUA,F,BDiffErr,CreateNewFigure=false);
    hold on 
    plot(F.x(I)/CtrlVar.PlotXYscale,F.y(I)/CtrlVar.PlotXYscale,".k",MarkerSize=3,DisplayName="Meas. locations")
    % hold on ; PlotGroundingLines(); PlotCalvingFronts();
    title("$(B_{\mathrm{Retrieved}}-B_{\mathrm{Meas}})/B_{\mathrm{Error}}$",Interpreter="latex")
    subtitle("")
    lg=legend;
    lg.String{1}="$(B_{\mathrm{Retrieved}}-B_{\mathrm{Meas}})/B_{\mathrm{Error}}$";
    lg.Interpreter="latex";


    ax5=nexttile;
    UaPlots(CtrlVar,MUA,F,Priors.B,CreateNewFigure=false);
    title("B prior") ;  subtitle("");

    ax6=nexttile;
    UaPlots(CtrlVar,MUA,F,InvFinalValues.B-Priors.B,CreateNewFigure=false);
    title("Retrieved B -  Prior B ")   ; subtitle("");


    % feels a bit clumsy way of ensuring that each tile as its own colorbar, but I can't think of a simpler approach
    set(figB,CurrentAxes=ax1) ;
    cl=clim;
    if min(cl) <0 && max(cl)> 0
        CM=cmocean('balanced',25,'pivot',0) ; colormap(ax1,CM);
    else
        CM=cmocean('balanced',25) ; colormap(ax1,CM);
    end

    set(figB,CurrentAxes=ax2) ;
    cl=clim;
    if min(cl) <0 && max(cl)> 0
        CM=cmocean('balanced',25,'pivot',0) ; colormap(ax2,CM);
    else
        CM=cmocean('balanced',25) ; colormap(ax2,CM);
    end

    set(figB,CurrentAxes=ax3) ;

    cl=clim;
    if min(cl) <0 && max(cl)> 0
        CM=cmocean('balanced',25,'pivot',0) ; colormap(ax3,CM);
    else
        CM=cmocean('balanced') ; colormap(ax3,CM);
    end

    set(figB,CurrentAxes=ax4) ;
    cl=clim;
    if min(cl) <0 && max(cl)> 0
        CM=cmocean('balanced',25,'pivot',0) ; colormap(ax4,CM);
    else
        set(figB,CurrentAxes=ax4) ;  CM=cmocean('balanced',25) ; colormap(ax4,CM);
    end

    set(figB,CurrentAxes=ax5) ;

    cl=clim;
    if min(cl) <0 && max(cl)> 0
        CM=cmocean('balanced',25,'pivot',0) ; colormap(ax5,CM);
    else
        CM=cmocean('balanced',25) ; colormap(ax5,CM);
    end

    set(figB,CurrentAxes=ax6) ;
    cl=clim;
    if min(cl) <0 && max(cl)> 0
        CM=cmocean('balanced',25,'pivot',0) ; colormap(ax6,CM);
    else
        CM=cmocean('balanced',25) ; colormap(ax6,CM);
    end


    % figB.Position=[200 200 1300 800];
    TB.TileSpacing="tight";
    TB.Padding="tight";
    %colormap(othercolor("Mdarkterrain",32))

end

%%
figh=FindOrCreateFigure("h retrieved"); clf(figh)

Th=tiledlayout("flow") ;

axh1=nexttile;
UaPlots(CtrlVar,MUA,F,F.h,CreateNewFigure=false);
title("Retrieved ice thickness") ; subtitle("");
CM=cmocean('-ice',15) ; colormap(CM);
%set(gca,'ColorScale','log')

axh2=nexttile;
UaPlots(CtrlVar,MUA,F,F.s-Priors.B,CreateNewFigure=false);
title("s-B prior") ; subtitle("");
CM=cmocean('-ice',15) ; colormap(CM);
%set(gca,'ColorScale','log')

CbarLink=linkprop([axh1 axh2],'CLim') ; assignin('base','CbarLink_clim',CbarLink)

Th.TileSpacing="tight";
Th.Padding="tight";

%%

end