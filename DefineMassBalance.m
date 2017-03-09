function [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)



as=zeros(MUA.Nnodes,1)+0.3;


%rhofw=1000;
%L=3.34e5;
%rho=917;
%cw=3974 ;
%Gt=5e-2;
%Gs=1.e3-3;


%uH=u0*tanh(Hc/Hc0);
%Tzd=T0*(b-B)/zref;
%ab=rho*cw*Gt*uH.*Tzd/rhofw/L;

switch UserVar.MisExperiment
    
    case 'ice0'
        
        % basal melt always zero
        ab=zeros(MUA.Nnodes,1);
        
    case 'ice1ra'
        
        % basal melt for t<=100, then zero
        if time <=100
            
            Hc0=75;
            Omega=0.2 ;
            z0=-100;
            ab=-Omega*tanh((b-B)/Hc0).* max(z0-b,0);
            
        else
            ab=zeros(MUA.Nnodes,1);
        end
        
    case 'ice1rr'
        
        % basal melt applied at all times
        Hc0=75;
        Omega=0.2 ;
        z0=-100;
        ab=-Omega*tanh((b-B)/Hc0).* max(z0-b,0);
        
    case 'ice2ra'
        
        % basal metl at 100 m/a for x>48km for the first 100 years, then zero
        if time<=100
            
            ab=zeros(MUA.Nnodes,1);
            I=MUA.coordinates(:,1)>480e3;
            ab(I)=-100;
            
        else
            ab=zeros(MUA.Nnodes,1);
        end
        
    case 'ice2rr'
        
        ab=zeros(MUA.Nnodes,1);
        I=MUA.coordinates(:,1)>480e3;
        ab(I)=-100;
        
    otherwise
        
        error('case not found')
        
end


end

