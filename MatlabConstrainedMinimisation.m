function [Cest,AGlenEst,VUA,JoptVector]=MatlabConstrainedMinimisation(CtrlVar,MUA,sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
        s,S,B,h,VUA,Xint,Yint,xint,yint,DTxy,TRIxy,DTint,TRIint,...
        AGlen,C,Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
        n,m,alpha,rho,rhow,g,GF,Itime,JoptVector)
    
    if CtrlVar.InfoLevelAdjoint>=10;
        fprintf(' minimisation method MatlabConstrainedMinimisation \n ')
    end
    
    nIt=CtrlVar.MaxAdjointIterations;
    if CtrlVar.AdjointRestart==0;
        JoptVector=zeros(nIt+1,6)+NaN; iJ=0;
    else
        iJ=size(JoptVector,1)-1;
        JoptVector=[JoptVector;zeros(nIt,6)+NaN];
        if iJ==-1 ; JoptVector=zeros(nIt+1,6)+NaN; iJ=0; end
    end
    
    
    
    wMeasInt=Grid1toGrid2(DTxy,wMeas,Xint,Yint);
        
    iA=strfind(CtrlVar.AdjointGrad,'A'); iC=strfind(CtrlVar.AdjointGrad,'C');
    CtrlVar.isAgrad=~isempty(iA); CtrlVar.isCgrad=~isempty(iC);
    
    
    C0=C;  AGlen0=AGlen;
    
    upA=zeros(length(C),1)+CtrlVar.AGlenmax;  lowA=zeros(length(C),1)+CtrlVar.AGlenmin;
    upC=zeros(length(C),1)+CtrlVar.Cmax;      lowC=zeros(length(C),1)+CtrlVar.Cmin;

    Eps=CtrlVar.Cmin/1000;
    C0=kk_proj(C0,upC,lowC,Eps);
    Eps=CtrlVar.AGlenmin/1000;
    AGlen0=kk_proj(AGlen0,upA,lowA,Eps);

    AGlenEst=AGlen0 ; Cest=C0;
    
    
    
    %%
    if CtrlVar.isCgrad
        fprintf(' C gradient optimisation step \n ')
                         
        [u,v,lambdauv,K]=SSTREAM2dNR(CtrlVar,MUA,s,S,B,h,u,v,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime);
        
        func=@(C) CostFunctionValueAndGradient(CtrlVar,MUA,s,S,B,h,u,v,...
            AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,...
            LAdjoint,LAdjointrhs,lambdaAdjoint,DTxy,TRIxy,...
            uMeas,vMeas,wMeasInt,C_prior,AGlen_prior,Cd,CAGlen,CC,GF);
        
        Misfit0=func(C0);
        
        lb=zeros(length(C),1)+CtrlVar.Cmin;
        ub=zeros(length(C),1)+CtrlVar.Cmax;
        
        % active-set is not a large-scale algorithm
        
        
        fprintf(' Initial misfit is %-g \n ',Misfit0)
        
        TypicalX=mean(C0)+zeros(numel(C0),1); TolCon=CtrlVar.Cmin;
 
        
        % interior point
        options=optimset('Algorithm','interior-point','MaxIter',CtrlVar.MaxAdjointIterations+1,'Diagnostics','on',...
            'Display','iter-detailed',...
            'InitTrustRegionRadius',1e-10,'InitBarrierParam',1,'ObjectiveLimit',1e-20,'MaxProjCGIter',5000,...
            'Hessian',{'lbfgs',30},...  % interior point only options
            'GradObj','on','MaxFunEvals',5000,'ScaleProblem','obj-and-constr',...
            'TolFun',1e-5,'TolX',1e-10,'TolProjCG',1,...
            'TypicalX',TypicalX,'TolCon',TolCon,'PlotFcns',@optimplotfval,'OutputFcn',@FminConOutputFcn);

        
%         %SQP
%         options=optimset('Algorithm','sqp','ObjectiveLimit',Misfit0/10,...
%             'MaxIter',CtrlVar.MaxAdjointIterations,'Diagnostics','on',...
%             'Display','iter-detailed',...
%             'GradObj','on','MaxFunEvals',5,...
%             'ScaleProblem','obj-and-constr','TolFun',1e-6,'TolX',1e-6,...
%             'TypicalX',TypicalX,'TolCon',TolCon,'PlotFcns',@optimplotfval);
        
        
        
        [Cest,Misfit1,exitflag,output,lambda,grad,hessian]=fmincon(func,C0,[],[],[],[],lb,ub,[],options);
        %[Cest,Misfit1,exitflag,output]=fminunc(func,C0,options);
        %Misfit1=func(Cest);

        
        fprintf(' Final misfit is %-g \t and ratio between final and initial misfit is %-g \n ',Misfit1,Misfit1/Misfit0)
        figure  ; PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,Cest); title('Cest matlab') ; colorbar
        if Misfit1>Misfit0 ;
            fprintf(' Estimate did not reduce misfit. Estimate set to initial value \n ')
            Cest=C0;
        end
        
    end
    %%
    
     %%
    if CtrlVar.isAgrad
        fprintf(' A gradient optimisation step \n ')
        
      
        func=@(AGlen) CostFunctionValueAndGradient(CtrlVar,MUA,s,S,B,h,u,v,...
            AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,...
            LAdjoint,LAdjointrhs,lambdaAdjoint,DTxy,TRIxy,...
            uMeas,vMeas,wMeasInt,C_prior,AGlen_prior,Cd,CAGlen,CC,GF);
        
        Misfit0=func(AGlen0);
        fprintf(' Initial misfit is %-g \n ',Misfit0)
        
        lb=zeros(length(AGlen),1)+CtrlVar.AGlenmin;
        ub=zeros(length(AGlen),1)+CtrlVar.AGlenmax;
        
        % active-set is not a large-scale algorithm
        
        % Interior point
        TypicalX=mean(AGlen0)+zeros(numel(AGlen0),1); TolCon= CtrlVar.AGlenmin/1000;
        
        options=optimset('Algorithm','interior-point','MaxIter',CtrlVar.MaxAdjointIterations+1,'Diagnostics','on',...
            'Display','iter-detailed',...
            'InitTrustRegionRadius',1e-5,'InitBarrierParam',1e-1,'ObjectiveLimit',1e-20,'MaxProjCGIter',5000,'Hessian',{'lbfgs',30},...  % interior point only options
            'GradObj','on','MaxFunEvals',5000,...
            'ScaleProblem','obj-and-constr','TolFun',1e-20,'TolX',1e-20,...
            'TypicalX',TypicalX,'TolCon',TolCon,'PlotFcns',@optimplotfval);
                
        [AGlenEst,Misfit1,exitflag,output,lambda,MisfitGratient1,hessian]=fmincon(func,AGlen0,[],[],[],[],lb,ub,[],options);
        
    
        %[Misfit1,MisfitGradient]=func(AGlenEst);
        fprintf(' Final misfit is %-g \t and ratio between final and initial misfit is %-g \n ',Misfit1,Misfit1/Misfit0)
        figure  ; PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,AGlenEst); title('AGlenEst') ; colorbar
        
        if Misfit1>Misfit0 ;
            fprintf(' Estimate did not reduce misfit. Estimate set to initial value \n ')
            AGlenEst=AGlen0;
        end
    
    end
    %%
      
    
    function stop = FminConOutputFcn(x, optimValues, state)
        
        switch state
            case 'iter'
                iJ=iJ+1; 
                JoptVector(iJ,1) = optimValues.fval;
            otherwise
        end
        
        stop=false;
    end
    
end






