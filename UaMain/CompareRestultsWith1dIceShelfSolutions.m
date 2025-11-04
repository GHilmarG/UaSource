function CompareRestultsWith1dIceShelfSolutions(Experiment,coordinates,connectivity,u,v,s,b,S,B,time,dt,AGlen,C,n,m,rho,rhow,g,alpha,nip,DTxy,TRIxy,DTint,TRIint,CtrlVar);
    

    x=coordinates(:,1);
    
    x0=min(x) ; l=max(x); N=100; q=zeros(N,1);
    xplot=linspace(x0,l,N);
    for I=1:N
        %[UserVar,AGlen,n]=DefineAGlenDistribution(Experiment,coordinates,connectivity,s,b,h,S,B,rho,rhow,Itime,time,CtrlVar);
        q(I)=quadgk(@fun,x0,xplot(I));
    end
    h=s-b;
    
    uAna=(mean(rho)*(1-mean(rho)/rhow)*g*mean(h)/4)^n *q;
    
    figure 
    H0=1000;
    
    plot(xplot/H0,uAna,'b-');
    hold on ; plot(x/H0,u,'rx')
    legend('Analytical','Numerical')
    xlabel('x (km)') ; ylabel('u')
    

     
    function  AGlen=fun(x)
        
        AGlen=DefineAGlenDistribution([],[x(:) x(:)]);
        AGlen=AGlen';
        
    end
     
     
end

