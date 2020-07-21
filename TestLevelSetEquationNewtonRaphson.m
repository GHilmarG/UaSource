
%%
load TestLevelSetEquationNewtonRaphson.mat


ind=MUA.coordinates(:,1) < 93.1e3 ; LSF=zeros(MUA.Nnodes,1) ; LSF(ind)=1 ; LSF(~ind)=-1;
%[LSF,UserVar,RunInfo]=ReinitializeLevelSet(UserVar,RunInfo,CtrlVar,MUA,LSF)  ;
F0.LSF=LSF; F1.LSF=F0.LSF; 

%F0.ub=F0.ub*0+1d3 ; F0.vb=F0.vb*0+1d3 ; 
%F1.ub=F1.ub*0+1d3 ; F1.vb=F1.vb*0+1d3 ; 
% F1.c=F1.c*0 ; F0.c=F0.c*0 ; 
CtrlVar.dt=1; 
CtrlVar.InfoLevelNonLinIt=10 ; CtrlVar.doplots=1;

[UserVar,RunInfo,LSF1,lambda]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1);