

%% Check numerical lower limits on AGlen and C, and make sure the smooth clipping towards those limits take place over small but numerical reasonable range

function CtrlVar=SetSmoothFloorWidths(CtrlVar,F)


if ~isfield(CtrlVar,"CminWidthRelative")

    CtrlVar.CminWidthRelative     = 1e-10 ;
end

if ~isfield(CtrlVar,"AGlenminWidthRelative")
    CtrlVar.AGlenminWidthRelative = 1e-10 ;
end

% (20 Nov 2026) The widths are relative to the median, not the mean, of the fields. In inversions the estimates of A and C can span
% many orders of magnitude, and the mean is then dominated by a few extreme values: on restart, the widths (and with them the
% effective A and C at nodes close to the lower limits) changed, and the cost function at the same parameters with them (Pine
% Island, nx150: AGlenminWidth 1.6e-18 at the start of the inversion, 7.7e-16 on restart; J increased by 9%).
CtrlVar.CminWidth     = max(CtrlVar.Cmin,     CtrlVar.CminWidthRelative    *median(F.C,"omitnan")) ;
CtrlVar.AGlenminWidth = max(CtrlVar.AGlenmin, CtrlVar.AGlenminWidthRelative*median(F.AGlen,"omitnan")) ;

fprintf("\n For numerical reasons A and C must always be finite.\n ")
fprintf("     Cmin=%g      CminWidth=%g \n ",CtrlVar.Cmin,CtrlVar.CminWidth);
fprintf(" AGlenmin=%g  AGlenminWidth=%g \n ",CtrlVar.AGlenmin,CtrlVar.AGlenminWidth);

end


