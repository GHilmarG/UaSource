
function [UserVar,F,l,InvFinalValues,xAdjoint,yAdjoint,RunInfo,gammaAdjoint]=...
    QuasiNewtonInversion...
    (UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)



%[UserVar,Cest,AGlenEst,Info,ub,vb,ud,vd,l,xAdjoint,yAdjoint,gammaAdjoint]=...
%    QuasiNewtonInversion(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)
    
%[UserVar,Cest,AGlenEst,Info,ub,vb,ud,vd,l,xAdjoint,yAdjoint,gammaAdjoint]=...
    %UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info)

error('No longer used')
    
    
xAdjoint=[] ;yAdjoint=[];

iA=strfind(CtrlVar.AdjointGrad,'A');
iC=strfind(CtrlVar.AdjointGrad,'C');
isAgrad=~isempty(iA);
isCgrad=~isempty(iC);

F.AGlen=InvStartValues.AGlen;
F.C=InvStartValues.C;
F.n=InvStartValues.n;
F.m=InvStartValues.m;

if ~isCgrad
    fprintf('QuasiNewtonInversion currently only implemented for a C inversion\n')
    error('Change input')
end

nIt=CtrlVar.MaxAdjointIterations;

if  isempty(RunInfo.JoptVector)
    RunInfo.JoptVector=zeros(nIt,7)+NaN; iJ=0;
else
    iJ=size(RunInfo.JoptVector,1);
    RunInfo.JoptVector=[RunInfo.JoptVector;zeros(nIt,7)+NaN];
end

if ~isempty(CtrlVar.AdjointInitialSearchStepSize) && ~isnan(CtrlVar.AdjointInitialSearchStepSize) && isfinite(CtrlVar.AdjointInitialSearchStepSize)
    gammaAdjoint=CtrlVar.AdjointInitialSearchStepSize;
elseif ~isempty(InvStartValues.InitialSearchStepSize) && ~isnan(InvStartValues.InitialSearchStepSize) && isfinite(InvStartValues.InitialSearchStepSize)
    gammaAdjoint=InvStartValues.InitialSearchStepSize ;
else
    gammaAdjoint=1;
end

if gammaAdjoint==0 ; gammaAdjoint=1 ; end

dJdC=F.C*0;

F.C=-q*dJdC;
ObjFunc=@(q,F) CalcMisfitFunction(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Priors,Meas);


F=@(q,ub,vb,ud,vd) CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlenEst,C0-q*dJdC,n,m,alpha,rho,rhow,g,GF,Priors,Meas);


[J,ObjFuncTerms,F,l,dJdu,RunInfo]=ObjFunc(0,F);

% var her

%F=@(q,ub,vb,ud,vd) CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlenEst,Cest-q*dJdC,n,m,alpha,rho,rhow,g,GF,Priors,Meas);
%[J,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,l,dIdu,kv,rh,nlInfo]=F(0,ub,vb,ud,vd);



if iJ==0
    iJ=iJ+1;  Info.JoptVector(iJ,1)=J; Info.JoptVector(iJ,2)=Idata;
    Info.JoptVector(iJ,3)=IRegC; Info.JoptVector(iJ,4)=IRegAGlen;
    Info.JoptVector(iJ,5)=IBarrierC; Info.JoptVector(iJ,6)=IBarrierAGlen;
    Info.JoptVector(iJ,7)=gammaAdjoint;
    Info.AdjointGrad{iJ}=CtrlVar.AdjointGrad;
end

fprintf('\n +++++++++++ At start of inversion:  \t Cost=%-g \t Idata=%-g \t RegC=%-g \t BarrierC=%-g \t gamma=%-g \n \n',J,Idata,IRegC,IBarrierC,gammaAdjoint)

for iteration=1:nIt
    
    CtrlVar.Iteration=iteration;
    C0=Cest; J0=J; AGlen0=AGlenEst;
    
    switch CtrlVar.AdjointMinimisationMethod
        
        case{'FixPointEstimationOfSlipperiness','FixPointEstimationOfSlipperiness:HessianGuestimate'}
            
            dIdCdata=CtrlVar.MisfitMultiplier*Calc_FixPoint_deltaC(CtrlVar,MUA,C0,m,GF,ub,vb,Meas.us,Meas.vs);
            dIdCreg=Calc_dIregdC(CtrlVar,MUA,Priors.CovC,C0,Priors.C);
            dIdCbarrier=Calc_dIdCbarrier(CtrlVar,MUA,C0);
            dJdC=dIdCdata+dIdCreg+dIdCbarrier;
            
        case {'QuasiNewtonInversion','QuasiNewtonInversion:HessianGuesstimate'}
            
            [UserVar,dJdC,dJdAGlen,ub,vb,ud,vd,xAdjoint,yAdjoint,dIdCreg,dIdAGlenreg,dIdCdata,dIdAGlendata,dIdCbarrier,dIdAGlenbarrier,lambdaAdjoint]=...
                AdjointGradientNR2d(...
                UserVar,CtrlVar,MUA,BCs,BCsAdjoint,s,b,h,S,B,ub,vb,ud,vd,l,AGlen0,C0,n,m,alpha,rho,rhow,g,GF,Priors,Meas);
            
        otherwise
            error('what case')
            
    end
    
    if ~isempty(strfind(CtrlVar.AdjointMinimisationMethod,'HessianGuesstimate'))
        % Hessians
        %
        % if J(x)=I(x)+0.5 (x-xa)' inv(C) (x-xa) then
        %  J(x+dx)=I(x)+  Ix dx + 0.5 dx' Ixx dx + 0.5 (x+dx-xa)' inv(C) (x+dx-xa)
        %         =I(x)+ Ix dx + 0.5 dx' H dx + 0.5 (x-xa)' inv(C) (x-xa) + 0.5 dx' inv(C) (x-da) + 0.5 (x-xa)' inv(C) dx + 0.5 dx' inv(C) dx
        %         =I(x)+ 0.5 (x-xa)' inv(C) (x-xa)+ Ix dx + 0.5 dx' H dx + (x-xa)' inv(C) dx + 0.5 dx' inv(C) dx
        %  or dJ(x)= Ix dx + 0.5 dx' H dx + (x-xa)' inv(C) dx + 0.5 dx' inv(C) dx
        %
        %   setting dJ/dx=0 gives
        %  dJx= Ix + H dx + (x-xa)' inv(C)  + inv(C) dx=0
        % and therfore  ( H + inv(C)) dx= - (Ix + inv(C) (x-xa))
        %
        % if I(x)=0 then inv(C) dx= - inv(C) (x-xa) or dx=-(x-xa) which makes sense
        % and if inv(C)=0 then  H dx=-Ix  which is the usual Newton system expression
        %
        % If J(x)=I1(x) +I2(x) then simply
        %  (I1xx+I2xx) dx =-(I1x+I2x)
        %
        % this is the same as above with I1=I and I2=0.5 (x-xa)' inv(C) (x-xa)
        % and I2x=inv(C) (x-xa) and I2xx=inv(C)
        %
        
        if CtrlVar.InfoLevelAdjoint>=10
            fprintf('QuasiNewtonInversion with guesstimated Hessians\n')
        end
        ddIddCfp=Calc_FixPoint_ddIddC(CtrlVar,MUA,ub,vb,ud,vd,C0,GF);
        ddIregddC=Calc_ddIregddC(CtrlVar,MUA,Priors.CovC);
        ddIddCbarrier=Calc_ddIddCbarrier(CtrlVar,MUA,C0);
        %dJdC=(ddIddCfp+ddIregddC+ddIddCbarrier)\dJdC;
        dJdC=(dIdCdata+dIdCreg+dIdCbarrier)./diag(ddIddCfp+ddIregddC+ddIddCbarrier) ; % since I'm approximating all Hessians using diagonal matrices this is OK
    end
    
    %  Newton step (I'm obmitting the neg sign here because I then go in the neg direction)
    %dJdCTest=(ddIddCfp+ddIregddC+ddIddCbarrier)\(dIdCfp+dIregdC+dIdCbarrier);
    %dJdC=(dIdC+dIdCreg+dIdCbarrier)./diag(ddIddCfp+ddIregddC+ddIddCbarrier) ;
    %dJdC=dIdC+dIdCreg+dIdCbarrier;
    
    if any(isnan(dJdC))
        save TestSave
        error(' NaN in dIdC')
    end
    
    %C1=C0-gamma*dJdC; C1=kk_proj(C1,CtrlVar.Cmax,CtrlVar.Cmin);
    
    F=@(q,ub,vb,ud,vd) CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen0,C0-q*dJdC,n,m,alpha,rho,rhow,g,GF,Priors,Meas);
    [J1,Idata1,IRegC1,IRegAGlen1,IBarrierC1,IBarrierAGlen1,ub,vb,ud,vd,l,dJdu,kv,rh,nlInfo]=F(gammaAdjoint,ub,vb,ud,vd);
    
    if J1/J0<0.5
        fprintf(' Initial step accected with J0=%-g \t J1=%-g \t and J1/J0=%-g \n ',J0,J1,J1/J0)
        J=J1 ; Idata=Idata1 ; IRegC=IRegC1 ; IRegAGlen=IRegAGlen1;
        IBarrierC=IBarrierC1; IBarrierAGlen=IBarrierAGlen1;
        
    else
        %% backtracking/line search
        % dI/dgamma =d I(J0+gamma gradf)= gradf'*gradf
        bb=gammaAdjoint; fa=J0 ; fb=J1 ; 
        slope0=[]; 
        %slope0=-dJdC'*dJdC;
        
        F=@(q,ub,vb,ud,vd) CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlenEst,C0-q*dJdC,n,m,alpha,rho,rhow,g,GF,Priors,Meas);
        nOut=10; listInF=[1:4] ;  listOutF=[7:10];
        
        temp=CtrlVar.InfoLevelBackTrack;
        CtrlVar.InfoLevelBackTrack=CtrlVar.InfoLevelAdjoint;
        CtrlVar.BacktrackingGammaMin=CtrlVar.BacktrackingGammaMinAdjoint;
        [gammaAdjoint,fgamma,BackTrackingInfoVector,ArgOut{1:nOut-1}]=BackTracking(slope0,bb,fa,fb,F,CtrlVar,nOut,listInF,listOutF,ub,vb,ud,vd);
        CtrlVar.InfoLevelBackTrack=temp;
        
        ub=ArgOut{6} ; vb=ArgOut{7}; ud=ArgOut{8} ; vd=ArgOut{9};
        [J,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,l,dJdu,kv,rh,nlInfo]=F(gammaAdjoint,ub,vb,ud,vd);
        
        fprintf(' Backtracking returns gamma=%-g and fgamma=%-g, and fgamma/J=%-g \n',gammaAdjoint,fgamma,fgamma/J)
        
        Test=(fgamma-J)/J;
        
        if abs(Test)>0.1
            fprintf('fgamma%-g \t J=%-g \n ',fgamma,J)
            warning('QuasiNewton:Test','those should be equal')
        end
        
        % J and fgamma should be almost the save
        
    end
    
    Cest=C0-gammaAdjoint*dJdC; % update the estimate for C
    Cest=kk_proj(Cest,CtrlVar.Cmax,CtrlVar.Cmin);
    
    fprintf('\n +++++++++++ Iteration %-i \t Cost=%-g \t \t J/J0=%-g \t Idata=%-g \t IRegC=%-g \t IBarrierC=%-g \n \n ',...
        iteration+Info.InverseIterations,J,J/J0,Idata,IRegC,IBarrierC)
    
    iJ=iJ+1;  Info.JoptVector(iJ,1)=J; Info.JoptVector(iJ,2)=Idata;
    Info.JoptVector(iJ,3)=IRegC; Info.JoptVector(iJ,4)=IRegAGlen;
    Info.JoptVector(iJ,5)=IBarrierC; Info.JoptVector(iJ,6)=IBarrierAGlen;
    Info.JoptVector(iJ,7)=gammaAdjoint;
    Info.AdjointGrad{iJ}=CtrlVar.AdjointGrad;
    
    if  (iteration==nIt  || iteration==1 ) && CtrlVar.doplots==1 && CtrlVar.InfoLevelAdjoint > 10
        
        figure ;
        subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,dJdC);  title(sprintf('dJdC it=%-i',iteration)) ; axis equal tight ; colorbar
        subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,Cest); title(sprintf('C it=%-i',iteration)) ; axis equal tight ; colorbar
        subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,log10(Cest)); title(sprintf('log10(C) it=%-i',iteration)) ; axis equal tight ; colorbar
        subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,Cest-C0); title(sprintf('changes in C it=%-i',iteration)) ; axis equal tight ; colorbar
        
    end
    
    
    if gammaAdjoint==0 
        fprintf(' gamma returned equal to zero. line search has stagnated. breaking out \n')
        break
    end
    
    if CtrlVar.InfoLevelAdjoint > 1000
        
        filename=sprintf('%s-AdjointIteration%04i',CtrlVar.Experiment,iteration);
        fprintf('Saving current estimate of C in %s \n',filename)
        save(filename,'CtrlVar','MUA','Cest','dJdC')
    
        
    end
    
end

I=isnan(Info.JoptVector(:,1)) ; Info.JoptVector(I,:)=[];
Info.InverseIterations=Info.InverseIterations+iteration;


end
