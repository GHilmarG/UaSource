
function [UserVar,F,l,InvFinalValues,RunInfo]=...
    InvertForModelParameters(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)

InvFinalValues=InvStartValues;

if isempty(CtrlVar.Inverse.InitialLineSearchStepSize) ||  CtrlVar.Inverse.InitialLineSearchStepSize==0
    CtrlVar.Inverse.InitialLineSearchStepSize=InvStartValues.SearchStepSize;
end

%% Define inverse parameters and anonymous function returning objective function, directional derivative, and Hessian
%

switch upper(CtrlVar.Inverse.InvertFor)
    
    case 'C'
        p0=InvStartValues.C;
    case 'LOGC'
        p0=log10(InvStartValues.C);
    case {'A','LOGA'}
        p0=InvStartValues.AGlen;
end

[J0,Gradient,Hessian,JPGouts]=JGH(p0,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);

% The parameters passed in the anonymous function are those that exist at the time the anonymous function is created.
func=@(p) JGH(p,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);


%%

if CtrlVar.Inverse.TestAdjointGradient
    
    
    % Get the gradient using the adjoint method
    [J,Gradient,H]=func(p0);
    
    
    % calc brute force gradient
    iRange=1:numel(p0);
    dJ = CalcBruteForceGradient(func,p0,iRange);
    [Gradient(iRange) dJ(iRange)   Gradient(iRange)-dJ(iRange) Gradient(iRange)./dJ(iRange)]
    mean(Gradient-dJ)
    
    IFig0=figure('Name','Inversion','NumberTitle','off');
    
    subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,Gradient) ;
    hold on
    PlotMuaMesh(CtrlVar,MUA);
    title('Adjoint gradient')
    
    subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,dJ) ;
    hold on
    PlotMuaMesh(CtrlVar,MUA);
    title('Brute force gradient')
    
    
    subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,Gradient-dJ) ;
    hold on
    PlotMuaMesh(CtrlVar,MUA);
    title('Difference between adjoint and brute force derivatives')
    
    subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,Gradient./dJ) ;
    hold on
    PlotMuaMesh(CtrlVar,MUA);
    title('Ratio between adjoint and brute force derivatives')
    
    IFig0.Position=[948.43 41.571 1246.3 1115.4];
    
    
else
    
    
    %%
    
    switch CtrlVar.Inverse.MinimisationMethod
        
        case {'QuasiNewtonInversion','QuasiNewtonInversion:HessianGuesstimate'}
            
             [p,RunInfo]=QuasiNewtonInversion2(CtrlVar,func,p0,RunInfo,MUA);
             
%            [UserVar,F,l,InvFinalValues,xAdjoint,yAdjoint,RunInfo,gammaAdjoint]=...
%                 QuasiNewtonInversion...
%                 (UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);
            
            %[UserVar,Cest,AGlenEst,Info,ub,vb,ud,vd,l,xAdjoint,yAdjoint,gammaAdjoint]=QuasiNewtonInversion(...
            %UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info);
            
        case 'AdjointProjectedGradient'
            
            [UserVar,Cest,AGlenEst,Info,ub,vb,ud,vd,xAdjoint,yAdjoint,gammaAdjoint]=AdjointProjectedGradient(...
                UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info);
            
            
        case 'MatlabOptimizationToolbox'
            
            
            [p,RunInfo]=InversionUsingMatlabOptimizationToolbox3(CtrlVar,func,p0,RunInfo);
            
            
            
        otherwise
            error('what case? ')
    end
    
    
    switch upper(CtrlVar.Inverse.InvertFor)
        
        case 'C'
            InvFinalValues.C=p;
        case 'LOGC'
            InvFinalValues.C=10.^p;
        case {'A','AGLEN'}
            InvFinalValues.AGlen=p;
        case {'LOGA','LOGAGLEN'}
            InvFinalValues.AGlen=10.^p;
    end
    
    
    [InvFinalValues.C,iU,iL]=kk_proj(InvFinalValues.C,CtrlVar.Cmax,CtrlVar.Cmin);
    [InvFinalValues.AGlen,iU,iL]=kk_proj(InvFinalValues.AGlen,CtrlVar.AGlenmax,CtrlVar.AGlenmin);
    
    
end

[J,dIdp,Hessian,Outs,F]=JGH(p,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);
InvFinalValues.J=J;
InvFinalValues.I=Outs.I;
InvFinalValues.R=Outs.R;
InvFinalValues.dJdp=dIdp; 
InvFinalValues.dIdp=Outs.dIdp;
InvFinalValues.dRdp=Outs.dRdp; 
InvFinalValues.SearchStepSize=RunInfo.Inverse.StepSize(end); 
                              




end
