         
function [Cest,AGlenEst,JoptVector,ub,vb,ud,vd,lx,ly,gammaAdjoint]=AdjointNR2D(CtrlVar,MUA,JoptVector,...
                   s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,...
                   sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
                   AGlen_prior,CAGlen,C_prior,CC,b_prior,uError,vError,wError,...
                   Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
                   GF)
   
               
    disp(' AdjointNR2D')

    error('No longer used')
    
    
    switch lower(CtrlVar.MisfitFunction)
        case 'uvwdiscrete'
            Cd=sparse(1:2*MUA.Nnodes+nInt,1:2*MUA.Nnodes+nInt,[uError.^2;vError.^2;wError.^2],2*MUA.Nnodes+nInt,2*MUA.Nnodes+nInt);   % I=(B u-d)' M (B u -d )
        case {'uvdiscrete','uvintegral'}
            Cd=sparse(1:2*MUA.Nnodes,1:2*MUA.Nnodes,[uError.^2;vError.^2],2*MUA.Nnodes,2*MUA.Nnodes);   % I=(B u-d)' inv(M) (B u -d )
        otherwise
            error(' what case? ' )
    end
    
    
    
    %% Minimisation using steepest descent or conjugraded gradients. Gradients are projected to box constraints
    
   
    
    switch CtrlVar.AdjointMinimisationMethod
        
         
           case {'QuasiNewtonInversion','QuasiNewtonInversion:HessianGuesstimate','FixPointEstimationOfSlipperiness'}

               
               [Cest,AGlenEst,JoptVector,ub,vb,ud,vd,lx,ly,gammaAdjoint]=QuasiNewtonInversion(...
                   CtrlVar,MUA,JoptVector,...
                   s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,...
                   sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
                   AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
                   Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
                   GF);
               
               if CtrlVar.doplots
                   figure   % only plotting here for C
                   [It,~]=size(JoptVector);
                   semilogy(0:It-1,JoptVector(:,1),'-ro') ; hold on
                   semilogy(0:It-1,JoptVector(:,2),'-bx') ;
                   semilogy(0:It-1,JoptVector(:,3),'-g+') ;
                   semilogy(0:It-1,JoptVector(:,5),'-m^') ;
                   legend('Cost function','Data misfit','C Reg','C barrier')
                   xlabel('Iteration') ;
                   hold off
               end
               
               
        case 'AdjointProjectedGradient'

            [Cest,AGlenEst,JoptVector,ub,vb,ud,vd,lx,ly,gammaAdjoint]=AdjointProjectedGradient(...
                CtrlVar,MUA,JoptVector,...
                s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,...
                sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
                AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
                Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
                GF);
            
            

            %%
            if CtrlVar.doplots
                figure   % only plotting here for C
                [It,~]=size(JoptVector);
                semilogy(0:It-1,JoptVector(:,1),'-ro') ; hold on
                semilogy(0:It-1,JoptVector(:,2),'-bx') ;
                semilogy(0:It-1,JoptVector(:,3),'-g+') ;
                semilogy(0:It-1,JoptVector(:,5),'-m^') ;
                legend('Cost function','Data misfit','C Reg','C barrier')
                hold off
            end
            %%
            
            
            
        case 'MatlabConstrainedMinimisation'
            
            gammaAdjoint=NaN;
            
            [Cest,AGlenEst,ub,vb,ud,vd,JoptVector]=MatlabConstrainedMinimisation(CtrlVar,MUA,ub,vb,ud,vd,sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
                s,S,B,h,Xint,Yint,xint,yint,DTxy,TRIxy,DTint,TRIint,...
                AGlen,C,Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
                n,m,alpha,rho,rhow,g,GF,Itime,JoptVector);
            
            save('MatlabConstrainedMinimisation','JoptVector')
            
            if CtrlVar.doplots==1
                figure
                [It,~]=size(JoptVector);
                semilogy(0:It-1,JoptVector(:,1),'-ro') ; hold on
            end
            
            
            
        case 'ProjectedBFGS'
            
            [Cest,AGlenEst,JoptVector,ub,vb,ud,vd,lx,ly,gammaAdjoint]=ProjectedBFGS(...
                CtrlVar,MUA,JoptVector,...
                s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,...
                sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
                AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
                Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
                GF);

%             [Cest,AGlenEst,ub,vb,ud,vd,JoptVector]=ProjectedBFGS(sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
%                 s,S,B,h,ub,vb,ud,vd,coordinates,connectivity,Xint,Yint,xint,yint,Boundary,DTxy,TRIxy,DTint,TRIint,...
%                 nip,AGlen,C,Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
%                 n,m,alpha,rho,rhow,g,GF,CtrlVar,JoptVector);
            
            save('ProjectedBFGS','JoptVector')
            
            if CtrlVar.doplots==1
                figure
                [It,~]=size(JoptVector);
                semilogy(0:It-1,JoptVector(:,1),'-ro') ; hold on
                semilogy(0:It-1,JoptVector(:,2),'-bx') ;
                semilogy(0:It-1,JoptVector(:,3),'-g+') ;
                semilogy(0:It-1,JoptVector(:,4),'-m^') ;
                legend('Cost function','Data misfit','Reg','Barrier')
                hold off
            end
            
            
        otherwise
            error('what case? ')
    end
    
end
