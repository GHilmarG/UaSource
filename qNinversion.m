

function [C,AGlenEst,u,v,JoptVector]=quasiNewtonNinversion(sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
        s,S,B,h,u,v,coordinates,connectivity,Xint,Yint,xint,yint,Boundary,DTxy,TRIxy,DTint,TRIint,...
        nip,AGlen,C,Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
        n,m,alpha,rho,rhow,g,GF,CtrlVar,Itime,JoptVector)
    
    
    %save TestSave
    %error('sdfa')
    
    %%
    %load TestSave  ; C=C0;
    
    
    nIt=CtrlVar.MaxAdjointIterations;
    gammamin=1e-5;
    
    if CtrlVar.AdjointRestart==0;
        JoptVector=zeros(nIt+1,6)+NaN; iJ=0;
    else
        iJ=size(JoptVector,1)-1;
        JoptVector=[JoptVector;zeros(nIt,6)+NaN];
        if iJ==-1 ; JoptVector=zeros(nIt+1,6)+NaN; iJ=0; end
    end
    
    
    
    
    
    uEleMeas=Nodes2EleMean(connectivity,uMeas); vEleMeas=Nodes2EleMean(connectivity,vMeas); speedEleMeas=sqrt(uEleMeas.*uEleMeas+vEleMeas.*vEleMeas);
    wMeasInt=Grid1toGrid2(DTxy,wMeas,Xint,Yint);
    AGlen0=AGlen; AGlenEst=AGlen0 ; 
    gamma=1 ; 
   
    %[u,v,lambdauv,K]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen0,C0,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
    %fprintf(' solved forward problem using start values \n ')
    %[J0,Idata0,IRegC0,IRegAGlen0,~,IBarrierC0,IBarrierAGlen0]=MisfitFunction(u,v,[],uMeas,vMeas,wMeasInt,C0,C_prior,AGlen0,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
    
  
    func=@(C,AGlen,u,v) CalcMisfitFunction(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,...    
      [],uMeas,vMeas,wMeasInt,C_prior,AGlen_prior,Cd,CAGlen,CC,GF,CtrlVar);
    
    [J,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,u,v,dIdu]=func(C,AGlen,u,v);
    
    
    iJ=iJ+1;  JoptVector(iJ,1)=J; JoptVector(iJ,2)=Idata;
    JoptVector(iJ,3)=IRegC; JoptVector(iJ,4)=IRegAGlen;
    JoptVector(iJ,5)=IBarrierC; JoptVector(iJ,6)=IBarrierAGlen;
    
    fprintf('\n +++++++++++ Iteration %-i \t Cost=%-g \t Idata0=%-g \t IRegC=%-g \t IBarrierC=%-g \n \n',0,J,Idata,IRegC,IBarrierC)
    
    
    close all
    
    for Iteration=1:nIt
        
       C0=C; J0=J;
        % Fix point idea:  
        % tau= He (u/c)^(1/m) -> u= C (tau/He)^m 
        % ignore the dependency of tau on C
        % u_meas=Ctrue (tau/He)^m
        % u_num= C (tau/He)^m
        % Ctrue=C+(Ctrue-C)=C+deltaC=C+(u_meas-u_num) (tau/He)^{-m}
        % I want u=u_meas  ->  deltaC=(u_meas-u0) (tau/He)^{-m}
        %                             =(u_meas-u0) speed/C
        %
        % Compare to NR applied to cost function:
        % I= 1/2 (umeas-u)^2
        % dI/dC=dI/du du/dC
        %       =(umeas-u) du/dC
        %       =(umeas-u) (tau/He)^m   (if I assume that dtau/dC=0 then du/dC=(tau/He)^m)
        %
        % dI^2/dC^2 C = -du/dC  (tau/He)^m 
        %             = -(tau/He)^(2m)
        %   Delta C= -H\dI/dC= (tau/He)^(-m)  (same result as from fix-point idea)
        %
        
        %  tb=He .* C.^(-1/m).*speed.^(1/m);
        %  (tb/He)^m=  speed/C
        %
        % One problem with using Delta C = (umeas-u) (tau/He)^(-m) 
        %                                = (umeas-u) speed/C
        % is that tau u is not dependent on C where floating
        % we can think about this as if C is then infinitly large
        % and that the actual slipperiness is C/He
        % Hence: Delta C= He (umeas-u) speed/C
        %
        uEle=Nodes2EleMean(connectivity,u); vEle=Nodes2EleMean(connectivity,v); speedEle=sqrt(uEle.*uEle+vEle.*vEle);
       
        
        [tbx,tby,tb,beta2] = CalcBasalTraction(u,v,C,m,GF,connectivity,CtrlVar);
               
        %dIdC=-(speedEleMeas-speedEle)./(tb./GF.ele).^m; % this is not based on dIdC= dI/du  du/dC, this is an estimate of the (negative) Newton step -H\grad I
        
        dIdC=-(speedEleMeas-speedEle)./(speedEle./C); % this is not based on dIdC= dI/du  du/dC, this is an estimate of the (negative) Newton step -H\grad I
        % the problem with this expression is that u is not equal to c (tau/He)^m when He->0 
        dIdC=dIdC.*GF.ele;
        
        %% adjoint approach
        
         [fc,gc,Idata,IReg,IBarrier] = CostFunctionValueAndGradient(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,...
        AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,...
        LAdjoint,LAdjointrhs,lambdaAdjoint,DTxy,TRIxy,...
        uMeas,vMeas,wMeasInt,C_prior,AGlen_prior,Cd,CAGlen,CC,GF,CtrlVar);
    
        save TestSave
        error('fds')
         
        
        if any(isnan(dIdC))
            save TestSave
            error(' NaN in dIdC')
        end
        
        %I=GF.ele<0.1 ; dIdC(I)=0;  % don't change C where floating
        %dIdC=dIdC.*GF.ele;          % don't change C where floating
        
        if gamma<gammamin ; gamma=gammamin ; end 
        
        C1=C0-gamma*dIdC;   %
        C1=kk_proj(C1,CtrlVar.Cmax,CtrlVar.Cmin);

                
        %[u,v,lambdauv,K]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen0,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
        %fprintf(' solved forward problem using start values \n ')
        %[J1,Idata1,IRegC1,IRegAGlen1,~,IBarrierC1,IBarrierAGlen1]=MisfitFunction(u,v,[],uMeas,vMeas,wMeasInt,C,C_prior,AGlen0,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
        
        
        [J1,Idata1,IRegC1,IRegAGlen1,IBarrierC1,IBarrierAGlen1,u,v,dIdu]=func(C1,AGlen0,u,v);
        
        
                
        if J1/J0<0.5
            fprintf(' Initial step accected with J0=%-g \t J1=%-g \t and J1/J0=%-g \n ',J0,J1,J1/J0)
            C=C1;
        else
            %% backtracking/line search
            % d I/dgamma =d I(J0+gamma gradf)= gradf'*gradf
             b=gamma; fa=J0 ; fb=J1 ; gradf=dIdC ; slope0=[];
            
             F=@(q,u,v) func(C0-q*dIdC,AGlen,u,v); nOut=9; listInF=[1 2] ; listOutF=[7 8];
             

            [gamma,fgamma,InfoVector,ArgOut{1:nOut-1}]=BackTracking(slope0,b,fa,fb,F,CtrlVar,nOut,listInF,listOutF,u,v);
            
            u=ArgOut{listOutF(1)-1} ; v=ArgOut{listOutF(2)-1};
            
            fprintf(' Backtracking returns gamma=%-g and fgamma=%-g \n',gamma,fgamma)
            C=C0-gamma*dIdC;
            C=kk_proj(C,CtrlVar.Cmax,CtrlVar.Cmin);
            
        end
        
        %%
        [J,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,u,v,dIdu]=func(C,AGlen0,u,v);
        fprintf('\n +++++++++++ Iteration %-i \t Cost=%-g \t \t J/J0=%-g \t Idata0=%-g \t IRegC=%-g \t IBarrierC=%-g \n \n ',Iteration,J,J/J0,Idata,IRegC,IBarrierC)
        
        iJ=iJ+1;  JoptVector(iJ,1)=J; JoptVector(iJ,2)=Idata;
        JoptVector(iJ,3)=IRegC; JoptVector(iJ,4)=IRegAGlen;
        JoptVector(iJ,5)=IBarrierC; JoptVector(iJ,6)=IBarrierAGlen;
        
        if  Iteration==nIt  || Iteration==1 ;
            
            figure ; PlotElementBasedQuantities(coordinates,connectivity,dIdC);  title(sprintf('dIdC it=%-i',Iteration)) ; axis equal tight ; colorbar
             
            figure ; PlotElementBasedQuantities(coordinates,connectivity,tb);  title(sprintf('tb it=%-i',Iteration)) ; axis equal tight ; colorbar
            
            
            figure ; PlotElementBasedQuantities(coordinates,connectivity,speedEle); title(sprintf('speed it=%-i',Iteration)) ; axis equal tight ; colorbar
            figure ; PlotElementBasedQuantities(coordinates,connectivity,log10(speedEle)); title(sprintf('log10(speed) it=%-i',Iteration)) ; axis equal tight ; colorbar
            figure ; PlotElementBasedQuantities(coordinates,connectivity,speedEleMeas); title(sprintf('measured it=%-i',Iteration)) ; axis equal tight ; colorbar
            
            
            
            figure ; PlotElementBasedQuantities(coordinates,connectivity,C); title(sprintf('C it=%-i',Iteration)) ; axis equal tight ; colorbar
            figure ; PlotElementBasedQuantities(coordinates,connectivity,log10(C)); title(sprintf('log10(C) it=%-i',Iteration)) ; axis equal tight ; colorbar
            figure ; PlotElementBasedQuantities(coordinates,connectivity,C-C0); title(sprintf('changes in C it=%-i',Iteration)) ; axis equal tight ; colorbar
            
        end
        
        if gamma==0 ;
            fprintf(' gamma returned equal to zero. line search has stagnated. breaking out \n')
            break
        end
    end
    
    
end
