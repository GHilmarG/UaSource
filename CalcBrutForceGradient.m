function [dIdC,dIdAGlen]=...
        CalcBrutForceGradient(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen0,C0,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar,...
        uMeas,vMeas,wMeasInt,C_prior,AGlen_prior,Cd,CAGlen,CC,GF)
    
    
    iA=strfind(CtrlVar.AdjointGrad,'A'); iC=strfind(CtrlVar.AdjointGrad,'C'); isAgrad=~isempty(iA); isCgrad=~isempty(iC);
    
    Nele=length(connectivity) ; C0=zeros(Nele,1)+mean(C0) ; % create element based slipperiness
    
    [u0,v0,lambda0]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen0,C0,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
    [J0,Idata0,IRegC0,IRegAGlen0]=MisfitFunction(u0,v0,[],uMeas,vMeas,wMeasInt,C0,C_prior,AGlen0,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
    
    
    N=numel(C0);   dIdC=zeros(N,1) ; dIdAGlen=AGlen0*0;
    
    
    delta=1e-8*norm(C0);
    
    if isCgrad
        parfor I=1:N
            fprintf(' %i ',round(100*I/N))
            unit=zeros(N,1) ; unit(I)=1;
            Cpert=C0+delta*unit;
            
            [uPert,vPert,lambda]=SSTREAM2dNR(s,S,B,h,u0,v0,coordinates,connectivity,Boundary,nip,AGlen0,Cpert,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
            [J,Idata,IRegC,IRegAGlen]=MisfitFunction(uPert,vPert,[],uMeas,vMeas,[],Cpert,C_prior,AGlen0,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
            dIdC(I)=(Idata-Idata0)/delta;
            
        end
    end
    
    
    save TestSave dIdC 
end

