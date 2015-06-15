         
function  [InvFinalValues,ub,vb,ud,vd,ubvbLambda,udvdLambda,lx,ly,Info]=...
            InvertForModelParameters(...
            Experiment,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,ubvbLambda,udvdLambda,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info)
 
    InvFinalValues=InvStartValues;
     
    switch CtrlVar.AdjointMinimisationMethod
 
           case {'QuasiNewtonInversion','QuasiNewtonInversion:HessianGuesstimate','FixPointEstimationOfSlipperiness'}
               
               [Cest,AGlenEst,Info,ub,vb,ud,vd,ubvbLambda,udvdLambda,lx,ly,gammaAdjoint]=QuasiNewtonInversion(...
                   Experiment,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,ubvbLambda,udvdLambda,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info);
               
        case 'AdjointProjectedGradient'

            [Cest,AGlenEst,Info,ub,vb,ud,vd,lx,ly,gammaAdjoint]=AdjointProjectedGradient(...
                Experiment,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,ubvbLambda,udvdLambda,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info);
            
        case 'MatlabConstrainedMinimisation'
            
            error('InvertForModelParameters:OptionBroken','MatlabConstrainedMinimisation broken')
            gammaAdjoint=NaN;
            
            [Cest,AGlenEst,ub,vb,ud,vd,JoptVector]=MatlabConstrainedMinimisation(CtrlVar,MUA,ub,vb,ud,vd,sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
                s,S,B,h,Xint,Yint,xint,yint,DTxy,TRIxy,DTint,TRIint,...
                Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
                alpha,rho,rhow,g,GF,Itime,JoptVector);
            
            save('MatlabConstrainedMinimisation','JoptVector')
            
            if CtrlVar.doplots==1
                figure
                [It,~]=size(JoptVector);
                semilogy(0:It-1,JoptVector(:,1),'-ro') ; hold on
            end
            
            
            
        case 'ProjectedBFGS'
            
            error('InvertForModelParameters:ProjectedBFGS','ProjectedBFGS broken')
            [Cest,AGlenEst,JoptVector,ub,vb,ud,vd,lx,ly,gammaAdjoint]=ProjectedBFGS(...
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
