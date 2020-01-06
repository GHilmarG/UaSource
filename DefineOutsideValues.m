function OutsideValue=DefineOutsideValues(UserVar,CtrlVar,MUA,F,OutsideValue)
    
    
    OutsideValue.h=2*CtrlVar.ThickMin ;
    
    OutsideValue.s=mean(F.S)+OutsideValue.h*(1-mean(F.rho)/F.rhow);
    OutsideValue.b=OutsideValue.s-OutsideValue.h;
    
    OutsideValue.ub=0;
    OutsideValue.vb=0;
    
    OutsideValue.ud=0;
    OutsideValue.vd=0;
    
    OutsideValue.dubdt=0;
    OutsideValue.dvbdt=0;
    
    OutsideValue.duddt=0;
    OutsideValue.dvddt=0;
    
    OutsideValue.dhdt=0;
    
    
end
