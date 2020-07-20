function TestLevelSetEquation(DoPlots,FileForSaving,Nsteps,RunType,Restart)

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
%
%   TestLevelSetEquation([],[],2000)
%
%  
%   TestLevelSetEquation(false,"TestLevelSet-RTinf-FAB1-CFp2q4-RCS-",5,"-RTinf-FAB1-CFp2q4-RCS-",false)
%
%   TestLevelSetEquation(false,"TestLevelSet-RTinf-FAB0k1-CFp2q4-RCS-",5,"-RTinf-FAB1-CFp2q4-RCS-",false)
%
%
%   job=batch('TestLevelSetEquation',0,{false,"TestLevelSet-RTinf-FAB1-CFp2q4-RCS-",5,"-RTinf-FAB1-CFp2q4-RCS-",false})
%
%   job=batch('TestLevelSetEquation',0,{false,"TestLevelSet-RTinf-FAB0k1-CFp2q4-RCS-",50000,"-RTinf-FAB1-CFp2q4-RCS-",false})
%   job10=batch('TestLevelSetEquation',0,{false,"TestLevelSet-RTinf-FAB10-CFp2q4-RCS-",50000,"-RTinf-FAB1-CFp2q4-RCS-",false})
%%
if nargin==0 || isempty(DoPlots)
    DoPlots=true;
end

if nargin<2 || isempty(FileForSaving)
    
    FileForSaving="TestLevelSetEquation.mat";
end

if nargin<2 || isempty(Nsteps)
    Nsteps=200;
end

if nargin<4 || isempty(RunType)
    RunType="-RTinf-FAB0k1-CFp2q4-" ; 
    RunType="-RTinf-FAB0k01-CFp4q4-" ; 
    RunType="-RTinf-FAB0k001-CFp4q4-" ; 
    RunType="-RTinf-FAB0k0001-CFp4q4-" ; 
    RunType="-RTinf-FAB1-CFp4q4-" ; % p4 can reduce NR convergence
    RunType="-RTinf-FAB1-CFp2q4-" ; % p4 can reduce NR convergence
end

if nargin<5
    Restart=false;
end


if Restart
    load("RestartFileTestLevelSetEquation"+FileForSaving,'CtrlVar','UserVar','RunInfo','MUA','F0','F1','xc','xcVector','BCs','N0')
else
    load('RestartEx-Reinitialize-RT0k5-FAB0-CubicMF-LevelSetWithMeltFeedback-1dIceShelf-MBice0-SUPGtaus-Adapt1.mat',...
        'RunInfo','CtrlVarInRestartFile','UserVarInRestartFile','F','MUA','BCs')
    CtrlVar=CtrlVarInRestartFile;
    UserVar=UserVarInRestartFile;
    xcInitial=200e3;  % this is the initial calving front location
    F.LSF=xcInitial-MUA.coordinates(:,1) ;
    BCs.LSFFixedValue=BCs.LSFFixedValue*0+xcInitial;
    CtrlVar.dt=0.1; CtrlVar.time=0;
    clear LevelSetEquation
    F0=F ; F1=F;
    xcVector.time=[];
    xcVector.Value=[];
    N0=1;
    xc=xcInitial;
    CtrlVar.LSFslope=-1;
end
 
% dtcritical=CalcCFLdt2D(UserVar,RunInfo,CtrlVar,MUA,F) ;

% if CtrlVar.dt > dtcritical
%     error('fasd')
% end

% example: RunType="-RTinf-FAB0k1-CFp2q4-"
CtrlVar.LevelSetReinitializeTimeInterval=str2double(replace(extractBetween(RunType,"RT","-"),"k","."));
CtrlVar.LevelSetFAB=str2double(replace(extractBetween(RunType,"FAB","-"),"k","."));
CtrlVar.LevelSetFABCostFunction=extractBetween(RunType,"CF","-") ;

CtrlVar.LevelSetSolverMaxIterations=100;
CtrlVar.LSFDesiredWorkAndForceTolerances=[1000 1e-6];
CtrlVar.LSFDesiredWorkOrForceTolerances=[100 1e-6];
CtrlVar.LSFExitBackTrackingStepLength=1e-4;
CtrlVar.LSFAcceptableWorkAndForceTolerances=[inf 1e-3];
CtrlVar.LSFAcceptableWorkOrForceTolerances=[100 1e-3];






%%
% 1) Test using c=u and where d/dt should be zero
%    Results:
%
%   No reinit and FAB=0 is stable and shows d\phi/dt=0, looks very good
%   No reinit and FAB=1 is stable and shows d\phi/dt=0, looks very good
%
%   Reinit every 20yr and FAB=0, crashes after about 58 years, noticeable discontinuity at each reset
%   Reinit every 20yr and FAB=1, works but discontionous at reint times
%
%  2) c=c.*TH+u.*(1-TH) ; with c0=100km
%
%  3) c=c.*TH+u.*(1-TH) ; with c0=200km
%


[tAna,xcAna]=xcVersusTimeForAnalyticalOneDimentionalIceShelf();



xcVector.time =[xcVector.time(:) ; NaN(Nsteps,1)]; 
xcVector.Value=[xcVector.Value(:) ; NaN(Nsteps,1)]; 
for I=N0:(Nsteps+N0-1)
    
    % Define calving rate and velociy
    
    
    q=-2;
    k=86322275.9814533 ;
    [s,b,u]=AnalyticalOneDimentionalIceShelf(CtrlVar,MUA);
    h=s-b;
    c=k*h.^q;
    
    % J=abs(F0.LSF)>50e3 ; c(J)=u(J);
    TH = TopHatApprox(1/10e3,MUA.coordinates(:,1)-xc,50e3) ;
    c=c.*TH+u.*(1-TH) ;
    % c=u ; % TestIng
    
    F0.ub=u ; F0.vb=u*0;
    F1.ub=u ; F1.vb=u*0;
    
    F0.c=c;
    F1.c=c;
    
        
            
    
    [UserVar,RunInfo,F1.LSF,Mask]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1); 
    CtrlVar.time=CtrlVar.time+CtrlVar.dt ;
    F0=F1 ;
    
    [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F1.LSF,0); xc=mean(xc,'omitnan') ;
    xcVector.time(I)=CtrlVar.time;
    xcVector.Value(I)=xc;
    
    % find slope around xc
    ind=abs(MUA.coordinates(:,1)-xc)<25e3; 
    y=F1.LSF(ind) ; x=MUA.coordinates(ind,1) ;
    sol=[ones(numel(find(ind)),1) x(:) ]\y(:);
    fprintf('\n \n ++++++++++++++++++++++++++   Slope %f \n \n ',sol(2))
    if contains(RunType,"-RCS-")
        CtrlVar.LSFslope=sol(2) ;
    else
        CtrlVar.LSFslope=-1 ;
    end
    
    if DoPlots
        
        % plots and post-processing
        fig=FindOrCreateFigure("Level Set Profile"+RunType);
        clf(fig) ;
        yyaxis left
        plot(MUA.coordinates(:,1)/1000,F1.LSF/1000,'.k')
        ylim([-500 200])
        hold on
        
        
        plot([xc xc]/1000,[-500 200],'--k')
        ylabel('$\varphi$ (km)','interpreter','latex')
        yyaxis right
        hold on
        plot(MUA.coordinates(:,1)/1000,F1.c,'.r')
        plot(MUA.coordinates(:,1)/1000,F1.ub,'.c')
        legend('$\phi$','$\phi_0$','$c$','$u$','interpreter','latex')
        title(sprintf('time %f ',CtrlVar.time))
        ylabel('$u(x)$ and $c(x)$ (m/yr)','interpreter','latex')
        xlabel('$x$ (km)','interpreter','latex')
        
        % load TestLevelSetEquation.mat
        fig=FindOrCreateFigure("Zero level position with time"+RunType);
        plot(xcVector.time,xcVector.Value/1000,'o')
        hold on
        tt=xlim;
        hold on; plot(tAna,xcAna/1000,'r','LineWidth',2) ;
        xlim(tt)
        title(sprintf('time %f ',CtrlVar.time))
        ylabel('$x_c$ (km)','interpreter','latex')
        xlabel('$t$ (yr)','interpreter','latex')
        
        
    end
    
    if mod(I,1000) ==0
        fprintf('Saving files %s. \n',FileForSaving)
        save("RestartFileTestLevelSetEquation"+FileForSaving,'CtrlVar','UserVar','RunInfo','MUA','F0','F1','xc','xcVector','BCs','N0')
        save(FileForSaving,'xcVector','RunType','tAna','xcAna','CtrlVar')
    end
    
end
N0=I;

fprintf('Saving files %s',FileForSaving)
save("RestartFileTestLevelSetEquation"+FileForSaving,'CtrlVar','UserVar','RunInfo','MUA','F0','F1','xc','xcVector','BCs','N0')
save(FileForSaving,'xcVector','RunType','tAna','xcAna','CtrlVar')
fprintf('All done.\n')
%%
% 


