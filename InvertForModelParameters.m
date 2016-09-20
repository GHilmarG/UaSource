         
function  [UserVar,InvFinalValues,ub,vb,ud,vd,l,xAdjoint,yAdjoint,Info]=...
            InvertForModelParameters(...
            UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info)
        
    InvFinalValues=InvStartValues;
     
    switch CtrlVar.AdjointMinimisationMethod
 
           case {'QuasiNewtonInversion','QuasiNewtonInversion:HessianGuesstimate','FixPointEstimationOfSlipperiness'}
               
               [UserVar,Cest,AGlenEst,Info,ub,vb,ud,vd,l,xAdjoint,yAdjoint,gammaAdjoint]=QuasiNewtonInversion(...
                   UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info);
               
        case 'AdjointProjectedGradient'

            [UserVar,Cest,AGlenEst,Info,ub,vb,ud,vd,xAdjoint,yAdjoint,gammaAdjoint]=AdjointProjectedGradient(...
                UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info);
            
                    
        case {'MatlabOptimizationToolbox:fmincon','MatlabOptimizationToolbox:fminunc','MatlabOptimizationToolbox'}
            
            [UserVar,Cest,AGlenEst,Info,ub,vb,ud,vd,xAdjoint,yAdjoint,gammaAdjoint]=InversionUsingMatlabOptimizationToolbox(...
                UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info);
             
            
        case 'ProjectedBFGS'
            
            error('InvertForModelParameters:ProjectedBFGS','ProjectedBFGS broken')
            [Cest,AGlenEst,JoptVector,ub,vb,ud,vd,xAdjoint,yAdjoint,gammaAdjoint]=ProjectedBFGS(...
                CtrlVar,MUA,JoptVector,...
                s,b,h,S,B,ub,vb,ud,vd,rho,rhow,alpha,g,...
                sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
                AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
                Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
                GF);

%             [Cest,AGlenEst,ub,vb,ud,vd,JoptVector]=ProjectedBFGS(sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
%                 s,S,B,h,ub,vb,ud,vd,coordinates,connectivity,Xint,Yint,xint,yint,Boundary,DTxy,TRIxy,DTint,TRIint,...
%                 nip,AGlen,C,Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
%                 n,m,alpha,rho,rhow,g,GF,CtrlVar,JoptVector);
 
            
        otherwise
            error('what case? ')
    end
    
    InvFinalValues.AGlen=AGlenEst;
    InvFinalValues.C=Cest;
    InvFinalValues.InitialSearchStepSize=gammaAdjoint;
    
end
