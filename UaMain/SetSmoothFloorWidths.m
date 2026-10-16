

%% Check numerical lower limits on AGlen and C, and make sure the smooth clipping towards those limits take place over small but numerical reasonable range

function CtrlVar=SetSmoothFloorWidths(CtrlVar,F)


if ~isfield(CtrlVar,"CminWidthRelative")

    CtrlVar.CminWidthRelative     = 1e-10 ;
end

if ~isfield(CtrlVar,"AGlenminWidthRelative")
    CtrlVar.AGlenminWidthRelative = 1e-10 ;
end

CtrlVar.CminWidth     = max(CtrlVar.Cmin,     CtrlVar.CminWidthRelative    *mean(F.C)) ;
CtrlVar.AGlenminWidth = max(CtrlVar.AGlenmin, CtrlVar.AGlenminWidthRelative*mean(F.AGlen)) ;

fprintf("\n For numerical reasons A and C must always be finite.\n ")
fprintf("     Cmin=%g      CminWidth=%g \n ",CtrlVar.Cmin,CtrlVar.CminWidth);
fprintf(" AGlenmin=%g  AGlenminWidth=%g \n ",CtrlVar.AGlenmin,CtrlVar.AGlenminWidth);

end


