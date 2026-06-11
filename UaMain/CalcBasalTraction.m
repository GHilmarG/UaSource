





function [tbx,tby,tb,eta,HeInt] = CalcBasalTraction(CtrlVar,UserVar,MUA,F,options)



%%
%
%    [tbx,tby,tb,eta] = CalcBasalTraction(CtrlVar,UserVar,MUA,F,options)
%
% Calculates basal traction from basal velocity using the sliding law.
%
% Returns either nodal or integration point values depending on the values of the logical optional input variables
%
%   CalcNodalValues=[true|false] 
%   CalcIntegrationPointValues=[true|false]
%
% The calculation at integration points is fully consistent with the way basal traction is calculated internally in Ua.
%
% Also returns the effective viscosity at integration points, if CalcIntegrationPointValues=true; 
%
%
% Note: This can only be used to calculate basal traction when using the SSTREAM and the Hybrid flow approximation. This will
% not return correct results for the SSHEET approximation!
%
%
%%

arguments
    CtrlVar struct
    UserVar struct
    MUA     struct
    F       {mustBeA(F,{'struct','UaFields','numeric'})}

    options.CalcNodalValues  logical = false
    options.CalcIntegrationPointValues  logical = true
    options.CalvingFrontColor char = "b"
    options.GroundingLineColor char = "r"
    options.PlotResults logical = false;
    options.FigureTitle string = ""
    options.FigureName string = ""


end

eta=[];
HeInt=[]; 


if isempty(F.ub)
    tbx=[]; tby=[] ; tb=[];
    return
end


if options.CalcNodalValues

    [tbx,tby,tb] = CalcBasalTractionAtNodes(CtrlVar,UserVar,MUA,F) ;

    if options.PlotResults

        FindOrCreateFigure(" nodal traction "+options.FigureName) ;
        [cbar,~,Par]=QuiverColorGHG(F.x/CtrlVar.PlotXYscale,F.y/CtrlVar.PlotXYscale,tbx,tby) ;
        PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],color=options.GroundingLineColor);
        PlotCalvingFronts(CtrlVar,MUA,F,color=options.CalvingFrontColor);
        PlotMuaBoundary(CtrlVar,MUA,"k--");
        title(cbar,"(kPa)")
        title("Basal tractions at nodal points")


        UaPlots(CtrlVar,MUA,F,tb,FigureTitle=" magnitude of basal traction at nodal points ")
        title("magnitude of basal traction at nodal points"+options.FigureTitle)
        set(gca,'ColorScale','log') 

    end
end

% Note if both nodal and integration point values are calculated, only the integration point values are returned. 
if options.CalcIntegrationPointValues

    CtrlVar.uvhMatrixAssembly.ZeroFields=false;
    CtrlVar.uvhMatrixAssembly.Ronly=false;
    CtrlVar.OnlyCalcBasalDragAndEffectiveViscosity=true ;
    [tbx,tby,tb,eta,HeInt] = CalcBasalTractionAtIntegrationPoints(CtrlVar,UserVar,MUA,F,F) ;


    if options.PlotResults

        [F.xint,F.yint] = CalcIntegrationPointsCoordinates(MUA) ;
        fbt=FindOrCreateFigure(" integration points traction "+options.FigureName) ;  clf(fbt);

        if options.CalcNodalValues  % if nodal values were also calculated, make them comparable by using same scaling as before.
            Par.QuiverSameVelocityScalingsAsBefore=1;
        else
            Par=[];
        end
        cbar=QuiverColorGHG(F.xint/CtrlVar.PlotXYscale,F.yint/CtrlVar.PlotXYscale,tbx,tby,Par) ;
        hold on
        PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],color=options.GroundingLineColor);
        PlotCalvingFronts(CtrlVar,MUA,F,color=options.CalvingFrontColor);
        PlotMuaBoundary(CtrlVar,MUA,"k--");
        title(cbar,"(kPa)")
        title("Basal tractions at integration points"+options.FigureTitle)

        UaPlots(CtrlVar,MUA,F,tb,FigureTitle=" magnitude of basal traction at integration points")
        set(gca,'ColorScale','log') 
        title("magnitude of basal traction at integration points")

        cbar=UaPlots(CtrlVar,MUA,F,log10(eta),FigureTitle=" effective viscosity"+options.FigureName) ;
        title("Effective viscosity at integration points"+options.FigureTitle)
        title(cbar,"$\log_{10}(\eta)$",interpreter="latex")


    end
end



%%

end

%% local functions


function [tbx,tby,tb] = CalcBasalTractionAtNodes(CtrlVar,UserVar,MUA,F)

narginchk(4,4)

%%
%
% Calculates basal traction from basal velocity using the sliding law.
%
% Returns nodal values
%
% Note: There is a slight inconsistency with respect to how this is done
% internally in Ua in the sense that the floating mask is here evaluated at
% nodes, whereas internally this is done at integration points.
%
%
%%




hf=F.rhow*(F.S-F.B)./F.rho ;
He = HeavisideApprox(CtrlVar.kH,F.h-hf,CtrlVar.Hh0);  % 1
delta = DiracDelta(CtrlVar.kH,F.h-hf,CtrlVar.Hh0) ;

[tbx,tby] = ...
    BasalDrag(CtrlVar,MUA,He,delta,F.h,F.B,F.S-F.B,F.rho,F.rhow,F.ub,F.vb,F.C,F.m,F.uo,F.vo,F.Co,F.mo,F.ua,F.va,F.Ca,F.ma,F.q,F.g,F.muk,F.V0);

tb=sqrt(tbx.^2+tby.^2);




end



function [tbxInt,tbyInt,tb,etaInt,Heint] = CalcBasalTractionAtIntegrationPoints(CtrlVar,UserVar,MUA,F0,F1,l1,BCs1)

narginchk(5,7)

if nargin<7 || isempty(BCs1) || isnan(BCs1)
    BCs1=BoundaryConditions();
end

if nargin<6 || isempty(l1) || isnan(l)
    l1=UaLagrangeVariables;
end

%%
%
% Calculates basal traction from basal velocity using the sliding law at integration points
%
%
%%

RunInfo=[];
CtrlVar.MassBalanceGeometryFeedback=1;  % Here the implicit mass-balance feedback needs to be disabled, as otherwise there is a risk of recursive call to DefineMassBalance (thanks to Emily Hill for spotting this).
[~,~,~,~,tbxInt,tbyInt,etaInt,Heint]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;

tb=sqrt(tbxInt.^2+tbyInt.^2);


end