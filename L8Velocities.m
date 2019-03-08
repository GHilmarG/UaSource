function [uMeas,vMeas,Err]=L8Velocities(UserVar,CtrlVar,x,y)

narginchk(4,4)

persistent  FGu FGv FGerror

if isempty(FGu)
    
    locdir=pwd;

    cd(UserVar.MatlabInterpolants.Folder)
    fprintf('Loading L8-2015 velocity interpolants ')
    load('L8-2015-GriddedInterpolants-1000m.mat','FGu','FGv','FGerror')
    fprintf('done\n')
    cd(locdir)
    
end

uMeas=FGu(x,y);
vMeas=FGv(x,y);
Err=FGerror(x,y);

T=1e3;
I=(uMeas>T & vMeas > T )| Err>1e3;
uMeas(I)=0 ; vMeas(I)=0 ; Err(I)=1e9;
uMeas=double(uMeas) ; vMeas=double(vMeas); Err=double(Err);

% Thwaites ice shelf
%I=x>-1650e3 & x < -1500e3 & y>-540e3 & y<-410e3 & GF.node < 0.5 ;


end

