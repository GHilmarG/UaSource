





function TestSensitivityMatrixCalculations(CtrlVar,MUA,F,BCs,l)


%% Calculates sensitivity matrices and tests them against brute-force finite-differences calculations at selected nodes
%
%
% Test show excellent agreement (apart from not respecting BCs.)
%
% see also: dFuvdA.m, dFuvdC.m , duvdCFunc.m
%
%load("TestSave.mat","CtrlVar","MUA","F","BCs","l");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A

Node=1050;
[dudA,dvdA]=duvdAFunc(CtrlVar,MUA,F,BCs);

dudA=dudA(:,Node);
dvdA=dvdA(:,Node);


UaPlots(CtrlVar,MUA,F,dudA,FigureTitle="du/dA")
UaPlots(CtrlVar,MUA,F,dvdA,FigureTitle="dv/dA")
UaPlots(CtrlVar,MUA,F,[dudA dvdA],FigureTitle="A velocity response")

% test for one particular node and compare to finite-differences

UserVar=[]; RunInfo=UaRunInfo;
[UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);

% second-order centered finite differences
dA=F.AGlen(Node)*0.00001;
A0=F.AGlen;

Ap=A0;
Ap(Node)=Ap(Node)+dA;
F.AGlen=Ap;
[UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l); up=F.ub; vp=F.vb;
F.AGlen=A0;

Am=A0;
Am(Node)=Am(Node)-dA;
F.AGlen=Am;
[UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l); um=F.ub; vm=F.vb;
F.AGlen=A0;

dudApert=(up-um)/(2*dA) ;  dvdApert=(vp-vm)/(2*dA) ;

UaPlots(CtrlVar,MUA,F,dudApert,FigureTitle="du/dA finite differences")
UaPlots(CtrlVar,MUA,F,dvdApert,FigureTitle="dv/dA finite differences")
UaPlots(CtrlVar,MUA,F,[dudApert dvdApert],FigureTitle="finite differences A velocity response")
% plot comparison

figAgrad=FindOrCreateFigure("du/dA gradient test") ;  clf(figAgrad)
plot(dudA,dudApert,"or") ;
hold on
axis equal
AX=axis;
plot([min(dudA) max(dudA)],[min(dudA) max(dudA)],"--k") ;
axis equal tight ;
xlabel(" $du/dA$",Interpreter="latex")  ;
ylabel("Finite difference $du/dA$",Interpreter="latex")
ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
axis on ; axis equal tight ; box off
title("Comparision betweenadjoint and finite-differences gradient calculations")
set(gcf,'Color','white')

figAgrad=FindOrCreateFigure("dv/dA gradient test") ;  clf(figAgrad)
plot(dvdA,dvdApert,"or") ;
hold on
axis equal
plot([min(dvdA) max(dvdA)],[min(dvdA) max(dvdA)],"--k") ;
axis equal tight ;
xlabel(" $dv/dA$",Interpreter="latex")  ;
ylabel("Finite difference $dv/dA$",Interpreter="latex")
ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
axis on ; axis equal tight ; box off
title("Comparision betweenadjoint and finite-differences gradient calculations")
set(gcf,'Color','white')


FindOrCreateFigure("duA test histogram")
duA=dudA-dudApert;  dvA=dvdA-dvdApert;
histogram(duA)
xlabel("Difference between exact and finite-differnece $du/dA$",Interpreter="latex")
[norm(duA)/norm(dudA) norm(dvA)/norm(dvdA)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% C

%

Node=500;


[dudC,dvdC]=duvdCFunc(CtrlVar,MUA,F,BCs);
dudC=dudC(:,Node);  dvdC=dvdC(:,Node);

UaPlots(CtrlVar,MUA,F,dudC,FigureTitle="du/dC")
UaPlots(CtrlVar,MUA,F,dvdC,FigureTitle="dv/dC")
UaPlots(CtrlVar,MUA,F,[dudC dvdC],FigureTitle="C velocity response")

% finite differences approach
UserVar=[]; RunInfo=UaRunInfo;
[UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);

% first-order forward finite differences
u0=F.ub; v0=F.vb;
dC=F.C(Node)*0.00001;
F.C(Node)=F.C(Node)+dC;
[UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
up=F.ub; vp=F.vb;

dudCpert=(up-u0)/dC ;  dvdCpert=(vp-v0)/dC ;

UaPlots(CtrlVar,MUA,F,dudCpert,FigureTitle="du/dC pert")
UaPlots(CtrlVar,MUA,F,dvdCpert,FigureTitle="dv/dC pert")
UaPlots(CtrlVar,MUA,F,[dudCpert dvdCpert],FigureTitle="finite differences C velocity response")

% plot comparison
figCgrad=FindOrCreateFigure("du/dC gradient test") ;  clf(figCgrad)

plot(dudC,dudCpert,"or") ;
hold on
axis equal
plot([min(dudC) max(dudC)],[min(dudC) max(dudC)],"--k") ;
axis equal tight ;
xlabel(" $du/dC$",Interpreter="latex")  ;
ylabel("Finite difference $du/dC$",Interpreter="latex")
ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
axis on ; axis equal tight ; box off
title("Comparision betweenadjoint and finite-differences gradient calculations")
set(gcf,'Color','white')


figCgrad=FindOrCreateFigure("dv/dC gradient test") ;  clf(figCgrad)
plot(dvdC,dvdCpert,"or") ;
hold on
axis equal
plot([min(dvdC) max(dvdC)],[min(dvdC) max(dvdC)],"--k") ;
axis equal tight ;
xlabel(" $dv/dC$",Interpreter="latex")  ;
ylabel("Finite difference $dv/dC$",Interpreter="latex")
ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
axis on ; axis equal tight ; box off
title("Comparision betweenadjoint and finite-differences gradient calculations")
set(gcf,'Color','white')

duC=dudC-dudCpert;  dvC=dvdC-dvdCpert;
[norm(duC)/norm(dudC) norm(dvC)/norm(dvdC)]

end
%%