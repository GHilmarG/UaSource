




function   [F,Fm1]=UpdateFtimeDerivatives(UserVar,RunInfo,CtrlVar,MUA,F,F0,BCs,l)


narginchk(8,8)

% (9 Oct 2026) Fm1 holds the rates of the previous time step, ie the rates of F0, together with the time step over which these
% rates were calculated (Fm1.dtRates). Both are needed for the second-order explicit estimate, see
% ExplicitEstimationUsingBackwardDifferences.m. Fm1 is a UaFields object.
Fm1=UaFields;
Fm1.dhdt=F0.dhdt ;
Fm1.dubdt=F0.dubdt ; Fm1.dvbdt=F0.dvbdt;
Fm1.duddt=F0.duddt ; Fm1.dvddt=F0.dvddt;
Fm1.dtRates=F0.dtRates ;

% (9 Oct 2026) F.dt is the time step. (CtrlVar.dt is only a copy of F.dt, kept for compatibility.)
if F.dt==0
    F.dhdt=[];
    F.dubdt=[]; F.dvbdt=[];
    F.dsdt=[] ; F.dbdt=[];
    F.dtRates=NaN;
    return

end


F.dhdt=(F.h-F0.h)/F.dt;
F.dsdt=(F.s-F0.s)/F.dt;
F.dbdt=(F.b-F0.b)/F.dt;

F.dubdt=(F.ub-F0.ub)/F.dt ;
F.dvbdt=(F.vb-F0.vb)/F.dt;

F.duddt=(F.ud-F0.ud)/F.dt ;
F.dvddt=(F.vd-F0.vd)/F.dt;

F.dtRates=F.dt ;   % (9 Oct 2026) the time step over which the above rates (backward differences) were calculated


fprintf("\n     UpdateFtimeDerivatives [max(abs(F.dubdt)) max(abs(F.dvbdt)) max(abs(F.dhdt)) ]=[%f %f %f]\n",max(abs(F.dubdt)),max(abs(F.dvbdt)),max(abs(F.dhdt)))



if CtrlVar.inUpdateFtimeDerivatives.SetAllTimeDerivativesToZero

    F.dubdt=F.dubdt*0;
    F.dvbdt=F.dvbdt*0;
    F.duddt=F.duddt*0;
    F.dvddt=F.dvddt*0;
    F.dhdt=F.dhdt*0;


    fprintf("CtrlVar.inUpdateFtimeDerivatives.SetAllTimeDerivativesToZero=%i \n",...
        CtrlVar.inUpdateFtimeDerivatives.SetAllTimeDerivativesToZero)
    fprintf("UpdateFtimeDerivatives: After modification: [max(abs(F.dubdt)) max(abs(F.dvbdt)) max(abs(F.dhdt)) ]=[%f %f %f]\n",...
        max(abs(F.dubdt)),max(abs(F.dvbdt)),max(abs(F.dhdt)))

else

    if CtrlVar.inUpdateFtimeDerivatives.SetTimeDerivativesDowstreamOfCalvingFrontsToZero

        if  CtrlVar.LevelSetMethod || ~isempty(F.LSF)

            if ~isempty(F.LSF)

                I=F.LSF< 0  ;
                F.dubdt(I)=0; F.dvbdt(I)=0;
                F.duddt(I)=0; F.dvddt(I)=0;
                F.dhdt(I)=0;

                fprintf("CtrlVar.inUpdateFtimeDerivatives.SetTimeDerivativesDowstreamOfCalvingFrontsToZero=%i \n",...
                    CtrlVar.inUpdateFtimeDerivatives.SetTimeDerivativesDowstreamOfCalvingFrontsToZero)
                fprintf("UpdateFtimeDerivatives: After modification: [max(abs(F.dubdt)) max(abs(F.dvbdt)) max(abs(F.dhdt)) ]=[%f %f %f]\n",...
                    max(abs(F.dubdt)),max(abs(F.dvbdt)),max(abs(F.dhdt)))

            end

        end

        if CtrlVar.inUpdateFtimeDerivatives.SetTimeDerivativesAtMinIceThickToZero

            I= (F.h <= 2*CtrlVar.ThickMin) | (F0.h <= 2*CtrlVar.ThickMin) ;
            F.dubdt(I)=0; F.dvbdt(I)=0;
            F.duddt(I)=0; F.dvddt(I)=0;
            F.dhdt(I)=0;

               fprintf("CtrlVar.inUpdateFtimeDerivatives.SetTimeDerivativesAtMinIceThickToZero=%i \n",...
                CtrlVar.inUpdateFtimeDerivatives.SetTimeDerivativesAtMinIceThickToZero)
            fprintf("UpdateFtimeDerivatives: After modification: [max(abs(F.dubdt)) max(abs(F.dvbdt)) max(abs(F.dhdt)) ]=[%f %f %f]\n",...
                max(abs(F.dubdt)),max(abs(F.dvbdt)),max(abs(F.dhdt)))
            

        end

    end


   

end

% if max(abs(F.dhdt)) >1e8 || max(abs(F.dubdt)) >1e8 ||  max(abs(F.dvbdt)) >1e8
%     fprintf("Check: [max(abs(F.dubdt)) max(abs(F.dvbdt))]=[%f %f]\n",max(abs(F.dubdt)),max(abs(F.dvbdt)))
%      I=isoutlier(F.dubdt,'median',ThresholdFactor=1000); F.dubdt(I)=0; F.dvbdt(I)=0; F.duddt(I)=0; F.dvddt(I)=0;
%     fprintf("                      After removing outliers: [max(abs(F.dubdt)) max(abs(F.dvbdt))]=[%f %f]\n",max(abs(F.dubdt)),max(abs(F.dvbdt)))
% end

%% comparison with other dh/dt estimates





end