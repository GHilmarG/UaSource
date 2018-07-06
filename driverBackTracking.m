

 func=@(x)   (x*0.11-1)^2+1  ;
 
 slope0=[] ;  b=[]  ; fa=[] ; fb=[]  ; 
 CtrlVar.BackTrackPlotResults=1 ;
 CtrlVar.BackTrackBeta=1000;
 [gamma,fgamma,InfoVector]=BackTracking(slope0,b,fa,fb,func,CtrlVar);
 
 
 
 