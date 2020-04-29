function [UserVar,RunInfo,LSF1,lambda]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCsLevelSet,F0,F1)
    %%
    %
    %
    %  %  df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f0)
    %
    
    persistent iCalls
    
    
    if ~CtrlVar.LevelSetMethod
        LSF1=F1.LSF;
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
    
    MLC=BCs2MLC(CtrlVar,MUA,BCsLevelSet);
    L=MLC.hL ; Lrhs=MLC.hRhs ;
    
    
    
    NRitMax=5;
    NRit=0 ;
    rWorkTolerance=1e-15;
    rForceTolerance=1e-15;
    while true
        
        NRit=NRit+1 ;
        
        [UserVar,K,R]=LevelSetEquationAssemblyNR(UserVar,CtrlVar,MUA,F0.LSF,F0.c,F0.ub,F0.vb,F1.LSF,F1.c,F1.ub,F1.vb);
        
        [dLSF,lambda]=solveKApe(K,L,R,Lrhs,[],[],CtrlVar);
        dLSF=full(dLSF);
        F1.LSF=F1.LSF+dLSF;
        
        
        % Newton decrement:  (R+L'lambda) dLSF
        %  R has the units area
        %  dLSF is dimentionless
        %
        % When calculating residuals make R dimentionless as well by dividing by area
        R=R/MUA.Area;
        if ~isempty(L)
            rWork=abs((R+L'*lambda)*dLSF) ;
            rForce=full((R+L'*lambda)'*(R+L'*lambda)) ;
        else
            rWork=abs(R'*dLSF);
            rForce=full(R'*R);
        end
        rDist=dLSF'*dLSF;
        
        
        %% more info
        fprintf('NR it %i  : \t rWork=%-15g \t rForce=%-15g \t rDist=%-15g  \n',NRit,rWork,rForce,rDist) ;
        
        FindOrCreateFigure('LSF1'); PlotMeshScalarVariable(CtrlVar,MUA,F1.LSF); title('LSF1')
        FindOrCreateFigure('LSF0');  PlotMeshScalarVariable(CtrlVar,MUA,F0.LSF); title('LSF0')
        FindOrCreateFigure('dLSF');  PlotMeshScalarVariable(CtrlVar,MUA,dLSF); title('dLSF1')
        prompt = 'Do you want more? Y/N [Y]: ';
%          str = input(prompt,'s');
%         if isempty(str)
%             str = 'Y';
%         end
        %%
        
        
        
        if NRit>NRitMax
            fprintf('LevelSetEquationNewtonRaphson: Maximum number of NR iterations (%i) reached. \n ',NRitMax)
            break
        end
        
        if (rWork < rWorkTolerance) || (rForce < rForceTolerance)
            fprintf('LevelSetEquationNewtonRaphson: NR iteration converged in %i iterations with rWork=%g and rForce=%g \n',NRit,rWork,rForce)
            break
        end
        
        
        
        
        
        
        
    end
    
    
    LSF1=F1.LSF ;
    
    
    
end

