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
c=f1*0 ;

x=MUA.coordinates(:,1) ;

% Define initial calving front

yc=linspace(min(MUA.coordinates(:,2)),max(MUA.coordinates(:,2)));
xc=yc*0+300e3 ;
X=[xc(:) yc(:) ] ; Y=[MUA.coordinates(:,1) MUA.coordinates(:,2)];
Dist= pdist2(X,Y,'euclidean','Smallest',1) ; Dist=Dist(:) ;
Outside=x > 300e3 ; % inside outside
Dist(Outside)=-Dist(Outside) ;
figure ; PlotMeshScalarVariable(CtrlVar,MUA,Dist) ;


GF.node=Dist ;
hold on
[xc,yc]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'r') ;



f1=Dist ;

CtrlVar.dt=1; CtrlVar.time=0;
Nsteps=100;
close
for I=1:Nsteps
    
    f0=f1;
    [UserVar,kv,rh]=LevelSetEquationAssembly(UserVar,CtrlVar,MUA,f0,c,u,v,kappa);
    
    L=[] ; Lrhs=[]  ;
    
    [f1,lambda]=solveKApe(kv,L,rh,Lrhs,[],[],CtrlVar);
    
    f1=full(f1);
    
    CtrlVar.time=CtrlVar.time+CtrlVar.dt ;
    
    figure(1) ;
    
    % Plots
    % plot(x/1000,f1,'r.')
    % PlotMeshScalarVariable(CtrlVar,MUA,f1);
    GF.node=f1 ;[xc,yc]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'r') ;
    
    xlim([300-10 300+200])
    title(sprintf('time %f ',CtrlVar.time))
    
    
    % re-initialize
    X=[xc(:) yc(:) ] ; Y=[MUA.coordinates(:,1) MUA.coordinates(:,2)];
    Dist= pdist2(X,Y,'euclidean','Smallest',1) ; Dist=Dist(:) ;
    Outside=
    Dist(Outside)=-Dist(Outside) ;
    f1=Dist ;
    
    
    
end

end

