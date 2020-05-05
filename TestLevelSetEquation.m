

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
x=MUA.coordinates(:,1) ;
y=MUA.coordinates(:,2) ;

% Define initial calving front as a (closed loop) boundary 
LargeDistance=1000e3 ;
xc0=300e3; 
CalvingFront=[xc0 min(y) ; xc0 max(y) ];
ymax=max(y) ; ymin=min(y) ; 
yc=linspace(ymin,ymax) ; 
xc=xc0+100e3*cos(2*2*pi*yc/(ymax-ymin));
CalvingFront=[xc(:) yc(:)] ; 


CalvingFrontClosure=[CalvingFront(end,:)+[0 10];   xc0 LargeDistance ;...
                     -LargeDistance LargeDistance ; -LargeDistance -LargeDistance ; ... 
                     xc0 -LargeDistance ; CalvingFront(1,:)+[0 -10]];  
CalvingFront=[CalvingFront ; CalvingFrontClosure ] ;

Npoints=1000 ; CalvingFront = interparc(Npoints,CalvingFront(:,1),CalvingFront(:,2),'linear'); % add some points


DistSigned=SignedDistance([x y],CalvingFront); 
figure ; PlotMeshScalarVariable(CtrlVar,MUA,DistSigned/1000) ;
hold on
plot(CalvingFront(:,1)/1000,CalvingFront(:,2)/1000,'-or')

colormap(othercolor('BuOr_12',1024)); ModifyColormap();
ModifyColormap(); 
%%
% TestIng
RunInfo=[] ;  BCsLevelSet=[]; 
F.c=zeros(MUA.Nnodes,1)-100000;
F.ub=u ; F.vb=v ;
F.LSF=DistSigned ;

MLC=BCs2MLC(CtrlVar,MUA,BCsLevelSet);
L=MLC.hL ; Lrhs=MLC.hRhs ;
CtrlVar.LevelSetMethod=true;



F1=F; F0=F; 
CtrlVar.LevelSetSolutionMethod="Newton-Raphson" ;
% CtrlVar.LevelSetSolutionMethod="Piccard" ;
fprintf('\n \n \n')
[UserVar,RunInfo,LSF1,lambda]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCsLevelSet,F0,F1); 
F.LSF ; 
%%
% re-initialize
%
% 1) mask
Mask=CalcMeshMask(CtrlVar,MUA,F.LSF,0); 
% 2) Distance from calving front

Dist=pdist2(MUA.coordinates(Mask.NodesOn,:),MUA.coordinates,'euclidean','Smallest',1) ;
Dist(Mask.NodesOut)=-Dist(Mask.NodesOut); 
figure ; PlotMeshScalarVariable(CtrlVar,MUA,Dist(:)) ;

% 3) Replace LSF with signed distance over In and Out nodes
F.LSF(Mask.NodesIn)=Dist(Mask.NodesIn) ;
F.LSF(Mask.NodesOut)=Dist(Mask.NodesOUt) ;

%%
CtrlVar.LevelSetSolutionMethod="Piccard" ;
[UserVar,RunInfo,LSF1Piccard,lambda]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCsLevelSet,F0); 

%%
F.LSF=DistSigned ;

CtrlVar.dt=1; CtrlVar.time=0;
Nsteps=100; InitializeInterval=10000;
close all

RunInfo=[] ;  BCsLevelSet=[]; 

for I=1:Nsteps
    
    
    
    F.c=zeros(MUA.Nnodes,1);
    F.ub=u ; F.vb=v ; 


    
    [UserVar,F.LSF,lambda]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCsLevelSet,F,F.LSF); 
    CtrlVar.time=CtrlVar.time+CtrlVar.dt ;
    
    
    
    if mod(I,InitializeInterval)==0
        % re-initialize
        
        [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F) ; 
        
        [ycMax,IycMax]=max(yc) ;
        [ycMin,IycMin]=min(yc) ;
        CalvingFrontClosure=[xc(IycMax) yc(IycMax)+10 ;   xc(IycMax) LargeDistance ; -1000e3 LargeDistance ; -1000e3 -1000e3 ; xc(IycMax) -1000e3 ; xc(IycMin) yc(IycMin)-10] ;
        CalvingFront=[xc(:) yc(:) ; CalvingFrontClosure ] ;
        Npoints=1000 ; CalvingFront = interparc(Npoints,CalvingFront(:,1),CalvingFront(:,2),'linear'); % add some points
        DistSigned=SignedDistance([x y],CalvingFront);
        F.LSF=DistSigned ;
        
    end
    %     figure(1) ;
    %
    %     % Plots
    %     % plot(x/1000,f1,'r.')
    %     hold off
    %     PlotMeshScalarVariable(CtrlVar,MUA,f1);
    %     hold on
    %     GF.node=f1 ;[xc,yc]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'r') ;
    %
    % xlim([300-10 300+200])
    
    
    figure(2) ; hold off
    plot(x/1000,F.LSF,'.')
    title(sprintf('time %f ',CtrlVar.time))
    
end



