function   [F,Fm1]=UpdateFtimeDerivatives(UserVar,RunInfo,CtrlVar,MUA,F,F0)

Fm1.dhdt=F0.dhdt ;
Fm1.dubdt=F0.dubdt ; Fm1.dvbdt=F0.dvbdt;
Fm1.duddt=F0.duddt ; Fm1.dvddt=F0.dvddt;

if CtrlVar.dt==0
    F.dhdt=[];
    F.dubdt=[]; F.dvbdt=[];
    F.dsdt=[] ; F.dbdt=[];
else
    F.dhdt=(F.h-F0.h)/CtrlVar.dt;
    F.dsdt=(F.s-F0.s)/CtrlVar.dt;
    F.dbdt=(F.b-F0.b)/CtrlVar.dt;
    
    F.dubdt=(F.ub-F0.ub)/CtrlVar.dt ; 
    F.dvbdt=(F.vb-F0.vb)/CtrlVar.dt;
    
    F.duddt=(F.ud-F0.ud)/CtrlVar.dt ; 
    F.dvddt=(F.vd-F0.vd)/CtrlVar.dt;
end



end