function [UserVar,f1,lambda]=LevelSetEquation(UserVar,CtrlVar,MUA,F,phi0,SF)

%%
% Solves the tracer conservation equation for the tracer c on the form:
% 
% $$\partial c/\partial t + d (u c)/dx + d (v c)/dy - \nabla \cdot (\kappa \nabla c) = a$$
% 
% The natural boundary condition is 
%
% $$\nabla c \cdot \hat{u} = 0 $$
%
% ie,  the free outflow condition 
%


%   < f | N + M >  


MLC=BCs2MLC(CtrlVar,MUA,BCsTracer);
L=MLC.hL ; Lrhs=MLC.hRhs ; lambda=Lrhs*0;

%%
Klear
load('RestartMismipPlus-ice0.mat','MUA','F','CtrlVarInRestartFile','BCs');
UserVar=[];
CtrlVar=CtrlVarInRestartFile ;

u=F.ub*0+1000 ; v=F.vb*0 ;
kappa=0 ;
f1=F.s*0+1;


x=MUA.coordinates(:,1) ;
y=MUA.coordinates(:,2) ;

% Define initial calving front as a (closed loop) boundary 
CalvingFront=[300e3 min(y) ; 300e3 max(y) ];

CalvingFrontClosure=[CalvingFront(end,:)+[0 10];   300e3 300e3 ; -300e3 300e3 ; -300e3 -300e3 ; 300e3 -300e3 ; CalvingFront(1,1)+[0 -10]];  
CalvingFront=[CalvingFront ; CalvingFrontClosure ] ;
Npoints=1000 ; CalvingFront = interparc(Npoints,CalvingFront(:,1),CalvingFront(:,2),'linear'); % add some points


DistSigned=SignedDistance([x y],CalvingFront); 
figure ; PlotMeshScalarVariable(CtrlVar,MUA,DistSigned/1000) ;
hold on
plot(CalvingFront(:,1)/1000,CalvingFront(:,2)/1000,'-or')

GF.node=DistSigned; 
hold on
[xc,yc]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'k') ;

%%

f1=DistSigned ;

CtrlVar.dt=1; CtrlVar.time=0;
Nsteps=10;
close all
figure
for I=1:Nsteps
    
    f0=f1;
    c=f0*0;
    [UserVar,kv,rh]=LevelSetEquationAssembly(UserVar,CtrlVar,MUA,f0,c,u,v,kappa);
    
    L=[] ; Lrhs=[]  ;
    
    [f1,lambda]=solveKApe(kv,L,rh,Lrhs,[],[],CtrlVar);
    
    f1=full(f1);
    
    CtrlVar.time=CtrlVar.time+CtrlVar.dt ;
    
    figure(1) ;
    
    % Plots
    % plot(x/1000,f1,'r.')
    hold off
    PlotMeshScalarVariable(CtrlVar,MUA,f1);
    hold on
    GF.node=f1 ;[xc,yc]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'r') ;
    
    % xlim([300-10 300+200])
    title(sprintf('time %f ',CtrlVar.time))
    
    
    % re-initialize
    
    [ycMax,IycMax]=max(yc) ; 
    [ycMin,IycMin]=min(yc) ;
    
    CalvingFrontClosure=[xc(IycMax) yc(IycMax)+10 ;   xc(IycMax) 300e3 ; -300e3 300e3 ; -300e3 -300e3 ; xc(IycMax) -300e3 ; xc(IycMin) yc(IycMin)-10] ;
    CalvingFront=[xc(:) yc(:) ; CalvingFrontClosure ] ;
    Npoints=1000 ; CalvingFront = interparc(Npoints,CalvingFront(:,1),CalvingFront(:,2),'linear'); % add some points
    DistSigned=SignedDistance([x y],CalvingFront);
    f1=DistSigned ;
    
   
    
end

end

