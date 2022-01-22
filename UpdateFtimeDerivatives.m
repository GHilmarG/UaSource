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


fprintf("\n     UpdateFtimeDerivatives [max(abs(F.dubdt)) max(abs(F.dvbdt))]=[%f %f]\n",max(abs(F.dubdt)),max(abs(F.dvbdt)))

if  CtrlVar.LevelSetMethod
    % I= (F.h <= 2*CtrlVar.LevelSetMinIceThickness) | (F0.h <= 2*CtrlVar.LevelSetMinIceThickness) ;
    if ~isempty(F.LSF)
        I=F.LSF< 0;
    end
else
    I= (F.h <= 2*CtrlVar.ThickMin) | (F0.h <= 2*CtrlVar.ThickMin) ;
end


F.dubdt(I)=0; F.dvbdt(I)=0;
F.duddt(I)=0; F.dvddt(I)=0;

fprintf("After removing values downstream of level set: [max(abs(F.dubdt)) max(abs(F.dvbdt))]=[%f %f]\n",max(abs(F.dubdt)),max(abs(F.dvbdt)))

I=isoutlier(F.dubdt,'median',ThresholdFactor=1000); F.dubdt(I)=0; F.dvbdt(I)=0; F.duddt(I)=0; F.dvddt(I)=0;
fprintf("                      After removing outliers: [max(abs(F.dubdt)) max(abs(F.dvbdt))]=[%f %f]\n",max(abs(F.dubdt)),max(abs(F.dvbdt)))

if max(abs(F.dubdt)) >1e8
    fprintf("Check: [max(abs(F.dubdt)) max(abs(F.dvbdt))]=[%f %f]\n",max(abs(F.dubdt)),max(abs(F.dvbdt)))
end

end