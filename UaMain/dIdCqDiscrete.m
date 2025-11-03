function [dIdC]=dIdCqDiscrete(u,v,lx,ly,S,B,h,C,m,rho,rhow,coordinates,connectivity,CtrlVar)
    
	     
	
	
	hf=rhow*(S-B)/rho;
    
	kH=CtrlVar.kH;
	Henod = HeavisideApprox(kH,h-hf);
	
	%C(C<CtrlVar.Cmin)=CtrlVar.Cmin;
    % There appears to be a potential problem when C for some reason is very small at one node
    % this can cause Ctemp to be very large and jumpy
    
    %C=EleAverageInterpolate(C,coordinates,connectivity);
    C(C<CtrlVar.CAdjointZero)=CtrlVar.CAdjointZero;
    
	Ctemp= (1/m)*Henod.*C.^(-1/m-1).*(sqrt(u.*u+v.*v+CtrlVar.SpeedZero^2)).^(1/m-1) ;
	dIdC=Ctemp.*(u.*lx+v.*ly);
	
%	AreaProNode=2e5*2e5/1089/1089;
%    dIdC=dIdC*AreaProNode;

end


