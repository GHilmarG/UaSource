function   [F,Fm1]=UpdateFtimeDerivatives(UserVar,RunInfo,CtrlVar,MUA,F,F0)

Fm1.dhdt=F0.dhdt ;
Fm1.dubdt=F0.dubdt ; Fm1.dvbdt=F0.dvbdt;
Fm1.duddt=F0.duddt ; Fm1.dvddt=F0.dvddt;

if CtrlVar.dt==0
    F.dhdt=[];
    F.dubdt=[]; F.dvbdt=[];
    F.dsdt=[] ; F.dbdt=[];
    return

end


F.dhdt=(F.h-F0.h)/CtrlVar.dt;
F.dsdt=(F.s-F0.s)/CtrlVar.dt;
F.dbdt=(F.b-F0.b)/CtrlVar.dt;

F.dubdt=(F.ub-F0.ub)/CtrlVar.dt ;
F.dvbdt=(F.vb-F0.vb)/CtrlVar.dt;

F.duddt=(F.ud-F0.ud)/CtrlVar.dt ;
F.dvddt=(F.vd-F0.vd)/CtrlVar.dt;


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

if max(abs(F.dhdt)) >1e8 || max(abs(F.dubdt)) >1e8 ||  max(abs(F.dvbdt)) >1e8
    fprintf("Check: [max(abs(F.dubdt)) max(abs(F.dvbdt))]=[%f %f]\n",max(abs(F.dubdt)),max(abs(F.dvbdt)))
     I=isoutlier(F.dubdt,'median',ThresholdFactor=1000); F.dubdt(I)=0; F.dvbdt(I)=0; F.duddt(I)=0; F.dvddt(I)=0;
    fprintf("                      After removing outliers: [max(abs(F.dubdt)) max(abs(F.dvbdt))]=[%f %f]\n",max(abs(F.dubdt)),max(abs(F.dvbdt)))
end

end