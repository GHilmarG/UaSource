function [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,F)


%%
% Defines mass balance along upper and lower ice surfaces.
%
%   [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)
%
%   [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,F);
%
%   [UserVar,as,ab,dasdh,dabdh]=DefineMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);
%
%   as        mass balance along upper surface 
%   ab        mass balance along lower ice surface
%   dasdh     upper surface mass balance gradient with respect to ice thickness
%   dabdh     lower surface mass balance gradient with respect to ice thickness
%  
% dasdh and dabdh only need to be specified if the mass-balance feedback option is
% being used. 
%
% In Úa the mass balance, as returned by this m-file, is multiplied internally by the local ice density. 
%
% The units of as and ab are water equivalent per time, i.e. usually
% as and ab will have the same units as velocity (something like m/yr or m/day).
%
% As in all other calls:
%
%  F.s       : is upper ice surface
%  F.b       : lower ice surface
%  F.B       : bedrock
%  F.S       : ocean surface
%  F.rhow    :  ocean density (scalar variable)
%  F.rho     :  ice density (nodal variable)
%  F.g       :  gravitational acceleration
%  F.x       : x nodal coordinates 
%  F.y       : y nodal coordinates 
%  F.GF      : The nodal grounded/floating mask (has other subfields)
%
%
%
% These fields need to be returned at the nodal coordinates. The nodal
% x and y coordinates are stored in MUA.coordinates, and also in F as F.x and F.y
%
% Examples:  
%
% To set upper surface mass balance to zero, and melt along the lower ice
% surface to 10 over all ice shelves:
%
%   as=zeros(MUA.Nnodes,1);
%   ab=-(1-F.GF.node)*10 
%
%
% To set upper surface mass balance as a function of local surface elevation and
% prescribe mass-balance feedback for the upper surface:
%
%   as=0.1*F.h+F.b;
%   dasdh=zeros(MUA.Nnodes,1)+0.1;
%   ab=F.s*0;
%   dabdh=zeros(MUA.Nnodes,1);
%
% 
%%


%rhofw=1000;
%L=3.34e5;
%rho=917;
%cw=3974 ;
%Gt=5e-2;
%Gs=1.e3-3;


%uH=u0*tanh(Hc/Hc0);
%Tzd=T0*(b-B)/zref;
%ab=rho*cw*Gt*uH.*Tzd/rhofw/L;


% zd is the ice draft, here it is b when afloat
% zb is the bedrock elevation, or B
% Hc=zd-db = b=B

as=0.3; 

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
            ab=-Omega*tanh((F.b-F.B)/Hc0).* max(z0-F.b,0);
            
            
            % when b>-100, for example b=-50 , then   z0-b=-100+50 =-50 < 0 
            
        else
            ab=zeros(MUA.Nnodes,1);
        end
        
    case 'ice1rr'
        
        % basal melt applied at all times
        Hc0=75;
        Omega=0.2 ;
        z0=-100;
        ab=-Omega*tanh((F.b-F.B)/Hc0).* max(z0-F.b,0);
        
    case 'ice2ra'
        
        % basal melt at 100 m/a for x>48km for the first 100 years, then zero
        if time<=100
            
            ab=zeros(MUA.Nnodes,1);
            I=F.x>480e3;
            ab(I)=-100;
            
        else
            ab=zeros(MUA.Nnodes,1);
        end
        
    case 'ice2rr'
        
        ab=zeros(MUA.Nnodes,1);
        I=F.x>480e3;
        ab(I)=-100;
        
    otherwise
        
        error('case not found')
        
end


end

