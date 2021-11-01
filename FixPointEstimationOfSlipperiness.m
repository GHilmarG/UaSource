function [Cest,AGlenEst,JoptVector,ub,vb,ud,vd,lx,ly]=FixPointEstimationOfSlipperiness(...=FixPointEstimationOfSlipperiness(...
    CtrlVar,MUA,JoptVector,...
    s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,...
    sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
    AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
    Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
    GF)

lx=[] ;ly=[];

iA=strfind(CtrlVar.AdjointGrad,'A'); iC=strfind(CtrlVar.AdjointGrad,'C'); isAgrad=~isempty(iA); isCgrad=~isempty(iC);

AGlenEst=AGlen ;
Cest=C;


if ~isCgrad
    fprintf('FixedPointEstimateOfSlipperiness only implemented for a C inversion\n')
    return
end

nIt=CtrlVar.MaxAdjointIterations;


if CtrlVar.AdjointRestart==0
    JoptVector=zeros(nIt+1,6)+NaN; iJ=0;
else
    iJ=size(JoptVector,1)-1;
    JoptVector=[JoptVector;zeros(nIt,6)+NaN];
    if iJ==-1 ; JoptVector=zeros(nIt+1,6)+NaN; iJ=0; end
end

uEleMeas=Nodes2EleMean(MUA.connectivity,uMeas); vEleMeas=Nodes2EleMean(MUA.connectivity,vMeas);
speedEleMeas=sqrt(uEleMeas.*uEleMeas+vEleMeas.*vEleMeas);

AGlen0=AGlen;

gamma=CtrlVar.AdjointInitialSearchStepSize;

[J,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,dIdu,K]=...
    CalcMisfitFunction(CtrlVar,MUA,s,b,h,S,B,ub,vb,ud,vd,AGlen0,n,Cest,m,rho,rhow,alpha,g,...
    sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
    AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
    Luv,Luvrhs,lambdauv,GF);


iJ=iJ+1;  JoptVector(iJ,1)=J; JoptVector(iJ,2)=Idata;
JoptVector(iJ,3)=IRegC; JoptVector(iJ,4)=IRegAGlen;
JoptVector(iJ,5)=IBarrierC; JoptVector(iJ,6)=IBarrierAGlen;

fprintf('\n +++++++++++ Iteration %-i \t Cost=%-g \t Idata0=%-g \t IRegC=%-g \t IBarrierC=%-g \n \n',0,J,Idata,IRegC,IBarrierC)


close all

for Iteration=1:nIt
    
    C0=Cest; J0=J;
    
    
    
    dIdC=Calc_FixPoint_deltaC(CtrlVar,MUA,C,m,GF,ub,vb,uMeas,vMeas);
    dIregdC=Calc_dIregdC(CtrlVar,CC,Cest,C_prior);
    dIdCbarrier=Calc_dIdCbarrier(CtrlVar,Cest);
    dJdC=dIdC+dIregdC+dIdCbarrier;
    
    if any(isnan(dJdC))
        save TestSave
        error(' NaN in dIdC')
    end
    
    %I=GF.ele<0.1 ; dIdC(I)=0;  % don't change C where floating
    %dIdC=dIdC.*GF.ele;          % don't change C where floating
    
    %if gamma<gammamin ; gamma=gammamin ; end
    
    C1=C0-gamma*dJdC;   %
    C1=kk_proj(C1,CtrlVar.Cmax,CtrlVar.Cmin);
    
    
    %[u,v,lambdauv,K]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen0,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
    %fprintf(' solved forward problem using start values \n ')
    %[J1,Idata1,IRegC1,IRegAGlen1,~,IBarrierC1,IBarrierAGlen1]=MisfitFunction(u,v,[],uMeas,vMeas,wMeasInt,C,C_prior,AGlen0,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
    
    
    %[J1,Idata1,IRegC1,IRegAGlen1,IBarrierC1,IBarrierAGlen1,VUA,dIdu]=func(C1,AGlen0,VUA);
    [J1,Idata1,IRegC1,IRegAGlen1,IBarrierC1,IBarrierAGlen1,ub,vb,ud,vd,dIdu,K]=...
        CalcMisfitFunction(CtrlVar,MUA,s,b,h,S,B,ub,vb,ud,vd,AGlen0,n,C1,m,rho,rhow,alpha,g,...
        sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
        AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
        Luv,Luvrhs,lambdauv,GF);
    
    
    
    if J1/J0<0.5
        fprintf(' Initial step accected with J0=%-g \t J1=%-g \t and J1/J0=%-g \n ',J0,J1,J1/J0)
        Cest=C1;
    else
        %% backtracking/line search
        % d I/dgamma =d I(J0+gamma gradf)= gradf'*gradf
        b=gamma; fa=J0 ; fb=J1  ; slope0=[];
        
        %F=@(q,VUA) func(C0-q*dJdC,AGlen,VUA); nOut=8; listInF=[1] ; listOutF=[7];
        F=@(q,ub,vb,ud,vd) CalcMisfitFunction(CtrlVar,MUA,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C0-q*dJdC,m,rho,rhow,alpha,g,...
            sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
            AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
            Luv,Luvrhs,lambdauv,GF);
        nOut=10;   % number of output arguments of func collected in BackTracking.m
        listInF=[1:4] ;
        listOutF=[7:10];  % listInF=1:10-7+1
        
        temp=CtrlVar.InfoLevelBackTrack;
        CtrlVar.InfoLevelBackTrack=CtrlVar.InfoLevelAdjoint;
        [gamma,fgamma,InfoVector,ArgOut{1:nOut-1}]=BackTracking(slope0,b,fa,fb,F,CtrlVar,nOut,listInF,listOutF,ub,vb,ud,vd);
        CtrlVar.InfoLevelBackTrack=temp;
        
        ud=ArgOut{6} ; vd=ArgOut{7}; ub=ArgOut{8} ; vb=ArgOut{9};
        %ub=ArgOut{listOutF(1)-1} ; vb=ArgOut{listOutF(2)-1};
        
        fprintf(' Backtracking returns gamma=%-g and fgamma=%-g \n',gamma,fgamma)
        Cest=C0-gamma*dJdC;
        Cest=kk_proj(Cest,CtrlVar.Cmax,CtrlVar.Cmin);
        
    end
    
    %%
    %[J,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,VUA,dIdu]=func(C,AGlen0,VUA);
    [J,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,dIdu,K]=...
        CalcMisfitFunction(CtrlVar,MUA,s,b,h,S,B,ub,vb,ud,vd,AGlen0,n,Cest,m,rho,rhow,alpha,g,...
        sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
        AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
        Luv,Luvrhs,lambdauv,GF);
    fprintf('\n +++++++++++ Iteration %-i \t Cost=%-g \t \t J/J0=%-g \t Idata0=%-g \t IRegC=%-g \t IBarrierC=%-g \n \n ',Iteration,J,J/J0,Idata,IRegC,IBarrierC)
    
    iJ=iJ+1;  JoptVector(iJ,1)=J; JoptVector(iJ,2)=Idata;
    JoptVector(iJ,3)=IRegC; JoptVector(iJ,4)=IRegAGlen;
    JoptVector(iJ,5)=IBarrierC; JoptVector(iJ,6)=IBarrierAGlen;
    
    if  Iteration==nIt  || Iteration==1 ;
        
        figure
        subplot(2,2,1) ; PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,dJdC);  title(sprintf('dJdC it=%-i',Iteration)) ; axis equal tight ; colorbar
        subplot(2,2,2) ; PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,speedEle); title(sprintf('speed it=%-i',Iteration)) ; axis equal tight ; colorbar
        subplot(2,2,3) ; PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,log10(speedEle)); title(sprintf('log10(speed) it=%-i',Iteration)) ; axis equal tight ; colorbar
        subplot(2,2,4) ; PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,speedEleMeas); title(sprintf('measured it=%-i',Iteration)) ; axis equal tight ; colorbar
        
        figure ;
        subplot(2,2,1) ; PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,dJdC);  title(sprintf('dJdC it=%-i',Iteration)) ; axis equal tight ; colorbar
        subplot(2,2,2) ; PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,Cest); title(sprintf('C it=%-i',Iteration)) ; axis equal tight ; colorbar
        subplot(2,2,3) ; PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,log10(Cest)); title(sprintf('log10(C) it=%-i',Iteration)) ; axis equal tight ; colorbar
        subplot(2,2,4) ; PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,Cest-C0); title(sprintf('changes in C it=%-i',Iteration)) ; axis equal tight ; colorbar
        
    end
    
    if gamma==0 ;
        fprintf(' gamma returned equal to zero. line search has stagnated. breaking out \n')
        break
    end
end

end
