function [UserVar,RunInfo,phi1,lambda]=LevelSetEquationPiccard(UserVar,RunInfo,CtrlVar,MUA,BCs,F0)
    %%
    %
    %
    %  %  df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f0)
    %
    % This was an initial rough approach, now only use LevelSetEquationNewtonRaphson.m
    %
    
    persistent iCalls
    
    
    if ~CtrlVar.LevelSetMethod
        phi1=F0.LSF; 
        lambda=[]; 
        return
    end

    
    if isempty(iCalls)
        iCalls=0 ;
    end
    iCalls=iCalls+1;
    
    
    % Define calving rate
    % F.c=zeros(MUA.Nnodes,1)-100e3 ;
    
    
    if ~isfield(CtrlVar,'LevelSetResetInterval') || isempty(CtrlVar.LevelSetResetInterval)
        CtrlVar.LevelSetResetInterval=10000;
    end
    
    
    MLC=BCs2MLC(CtrlVar,MUA,BCs);
    L=MLC.LSFL ; Lrhs=MLC.LSFRhs ;
    l=Lrhs*0;
  
    
    [UserVar,kv,rh]=LevelSetEquationAssembly(UserVar,CtrlVar,MUA,F0.LSF,F0.c,F0.ub,F0.vb);
    [phi1,lambda]=solveKApe(kv,L,rh,Lrhs,[],[],CtrlVar);
    phi1=full(phi1);
    
end

