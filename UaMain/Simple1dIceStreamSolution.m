function [uAna,xAna,GLpos,dudxGL,dl]=Simple1dIceStreamSolution(x,s,b,S,B,alpha,g,rho,rhow,A,C)
    
        
    
    if nargin==0
        n=1 ; m=1;
        h=1000; alpha=0.001; g=9.81/1000 ; rho=900; rhow=1000;
        SlipRatio=10; ud=0.1 ;
        tau=rho*g*h*sin(alpha);
        A=ud/(2*tau^n*h/(n+1));
        ub=SlipRatio*ud;
        C=ub/tau^m;
        L=100e3;
        x=0:1000:L;
        S=0; B0=0;
        B=B0-tan(alpha)*x;
        
        hf=rhow*(S-B)./rho ;
        gf =double(h>hf);  % 1 if grounded, 0 if afloat
        
        
        bfloat=S-rho.*h/rhow;
        
        b=gf.*B + (1-gf) .* bfloat ;
        s=gf.*(B+h) + (1-gf) .* (bfloat+h) ;
        
        
    end
    
    
    
    h=mean(s-b);
    A=mean(A);
    C=mean(C);
    rho=mean(rho);
    
    tau=rho*g*h*alpha;
    hf=rhow*(S-B)./rho;
    
    GLpos=(rho*h/rhow+B(1)-S(1))/tan(alpha);
    l=min([GLpos max(x)]);
    
    dl=S(1)-B(1)+l*tan(alpha);   % d at grounding line, or at end of domain if GL not in model domain
    dl=max([dl 0]);
   % gf =double(h>hf);  % 1 if grounded, 0 if afloat
    d=S-b; d(d<0)=0;

%     if ~any(gf==0)
%         l=max(x) ; dl=d(end);
%     else
%         [~,ind]=max(x.*gf);
%         %ind=find(diff(gf)==-1);
%         dl=d(ind) ;
%         l=x(ind);
%     end
    
    
    
    xAna=linspace(0,l,10000);
    kappa=sqrt(A/(2*h*C));
    
    n=1;
    K=A*(g*(rho*h^2 - rhow*dl^2)/4/h)^n;
    
    dudxGL=K;
    
%     %c2=(K*exp(-kappa*l)+C*tau*kappa)/(kappa*(exp(-2*kappa*l)-1));
%     c2=-(K*exp(-kappa*l)+C*tau*kappa)/(kappa*(exp(-2*kappa*l)+1));
%     c2=-(K+C*tau*kappa*exp(kappa*l))/(2*kappa*cosh(kappa*l));
%     c1=-c2-C*tau;
    
    c1=K/(2*kappa*cosh(kappa*l)); 
    c2=-c1;

   % uAna=c1*exp(kappa*xAna)+c2*exp(-kappa*xAna)+C*tau;

    uAna=C*tau+2*c1*sinh(kappa*xAna);
    uAna=C*tau+K*sinh(kappa*xAna)/(kappa*cosh(kappa*l));
    
    
    %dudx=gradient(uAna,xAna(2)-xAna(1));
    %K
    %dudx(end-3:end)
    
    
    if nargout==0
        figure ; plot(xAna/1000,s,'b') ; hold on ; plot(xAna/1000,d,'r') ; plot(xAna/1000,b,'g') ; legend('s','d','b')
        figure ; plot(xAna/1000,uAna)
    end
       
    
end


