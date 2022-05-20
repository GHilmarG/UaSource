function [rP,rL,rTP,rTL,rPTL]=CalcLSFconstfunctionTLPterms(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l,LSF)


F1.LSF=LSF ; 


% Testing residuals towards an automated criterion
gamma=0 ; 
L=[] ; Lrhs=[] ; dLSF=0 ; dl=0 ; 



CtrlVar.LSF.P=1 ; CtrlVar.LSF.T=0 ; CtrlVar.LSF.L=0 ; % Fixed point
[rP,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionLevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l,dLSF,dl,BCs) ; 

CtrlVar.LSF.P=1 ; CtrlVar.LSF.T=1 ; CtrlVar.LSF.L=0 ; % Pseudo time stepping
[rTP,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionLevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l,dLSF,dl,BCs) ; 

CtrlVar.LSF.P=1 ; CtrlVar.LSF.T=1 ; CtrlVar.LSF.L=1 ; % Full
[rPTL,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionLevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l,dLSF,dl,BCs) ;

CtrlVar.LSF.P=0 ; CtrlVar.LSF.T=0 ; CtrlVar.LSF.L=1 ; % advection term
[rL,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionLevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l,dLSF,dl,BCs) ;

CtrlVar.LSF.P=0 ; CtrlVar.LSF.T=1 ; CtrlVar.LSF.L=1 ; % advection term
[rTL,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionLevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l,dLSF,dl,BCs) ;

fprintf('      Diffusion term only, ie     P residual %g \n',rP)
fprintf('      Advection term only, ie     L residual %g \n',rL)
fprintf('  time and diffusion term, ie   T+P residual %g \n',rTP)
fprintf('  time and advection term, ie   T+L residual %g \n',rTL)
fprintf('            Full residual, ie T+L+P residual %g \n',rPTL)
fprintf('                                         P/L %g \n',rP/rL)




end
