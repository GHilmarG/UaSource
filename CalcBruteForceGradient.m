function [dJdC,dJdAGlen,dIdatadC,dIdatadAGlen,dIRegCdC,dIRegAGlendAGlen]=...
        CalcBruteForceGradient(CtrlVar,MUA,s,S,B,h,u,v,AGlen0,C0,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,...
        uMeas,vMeas,wMeasInt,C_prior,AGlen_prior,Cd,CAGlen,CC,GF)
    
    fprintf(' Calculating gradients using brute force method \n')
    
    iA=strfind(CtrlVar.AdjointGrad,'A'); iC=strfind(CtrlVar.AdjointGrad,'C'); isAgrad=~isempty(iA); isCgrad=~isempty(iC);
     
    C0=zeros(MUA.Nele,1)+mean(C0) ; % create element based slipperiness
    
    [u0,v0,lambda0]=SSTREAM2dNR(CtrlVar,MUA,s,S,B,h,u,v,AGlen0,C0,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g);
    [J0,Idata0,IRegC0,IRegAGlen0]=MisfitFunction(CtrlVar,MUA,u0,v0,[],uMeas,vMeas,wMeasInt,C0,C_prior,AGlen0,AGlen_prior,Cd,CAGlen,CC,GF);
    
    
    N=numel(C0);   dJdC=zeros(N,1) ; dJdAGlen=AGlen0*0; dIdatadC=zeros(N,1) ; dIdatadAGlen=AGlen0*0; dIRegCdC=zeros(N,1); dIRegAGlendAGlen=zeros(N,1);
    
    if isCgrad
        delta=1e-2*mean(C0);
        parfor I=1:N
            fprintf(' %i ',round(100*I/N))
            unit=zeros(N,1) ; unit(I)=1;
            Cpert=C0+delta*unit;
            
            [uPert,vPert,lambda]=SSTREAM2dNR(CtrlVar,MUA,s,S,B,h,u0,v0,AGlen0,Cpert,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g);
            [J,Idata,IRegC,IRegAGlen]=MisfitFunction(CtrlVar,MUA,uPert,vPert,[],uMeas,vMeas,[],Cpert,C_prior,AGlen0,AGlen_prior,Cd,CAGlen,CC,GF);
            dJdC(I)=(J-J0)/delta;
            dIdatadC(I)=(Idata-Idata0)/delta;
            dIRegCdC(I)=(IRegC-IRegC0)/delta;
            
        end
    end
    
    if isAgrad
        delta=1e-2*norm(AGlen0);
        parfor I=1:N
            fprintf(' %i ',round(100*I/N))
            unit=zeros(N,1) ; unit(I)=1;
            AGlenPert=AGlen0+delta*unit;
            
            [uPert,vPert,lambda]=SSTREAM2dNR(CtrlVar,MUA,s,S,B,h,u0,v0,AGlenPert,C0,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g);
            [J,Idata,IRegC,IRegAGlen]=MisfitFunction(CtrlVar,MUA,uPert,vPert,[],uMeas,vMeas,[],C0,C_prior,AGlenPert,AGlen_prior,Cd,CAGlen,CC,GF);
            dJdAGlen(I)=(J-J0)/delta;
            dIdatadAGlen(I)=(Idata-Idata0)/delta;
            dIRegAGlendAGlen(I)=(IRegAGlen-IRegAGlen0)/delta;
            
        end
    end
    
    save BrutForceGradientResults
    
end

