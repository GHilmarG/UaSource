function  [c,dcddphidx,dcddphidy]=DefineCalvingAtIntegrationPoints(UserVar,CtrlVar,dphidx,dphidy,F) 

% function  [c,dcDdphidx,dcDdphidy]=DefineCalvingAtIntegrationPoints(UserVar,CtrlVar,dfdx,dfdy,u,v,h,s,S,x,y)

%
%  phi is the level set function
%

%  dcddphidx is   dc/d (dphidx) = \frac{ d c]{d (dphidx)}, ie it is the derivative of c with respect to dphidx
%  where dphidx is in turn d(phi)/dx              
%
% cint=DefineCalvingAtIntegrationPoints(UserVar,CtrlVar,nx,ny,uint,vint)


isPlotStrainRates=false;

if isPlotStrainRates

    FigStrainRates=FindOrCreateFigure("strain rates and velocities") ; clf(FigStrainRates) ; 
    QuiverColorGHG(F.x/1000,F.y/1000,F.ub,F.vb) ;
    Scale=2;
    hold on ; PlotTensor(F.x/1000,F.y/1000,F.exx,F.exy,F.eyy,Scale) ; axis equal

end


narginchk(5,5)

switch UserVar.CalvingLaw.Type

    case "-AC-"

        CliffHeight=min((F.s-F.S),F.h) ;  % note: this does not account for density

        fI=3.2e-17*365.25 ; c=fI*CliffHeight.^(7.2) ;
        % Now set calving rate to zero for cliff less than 135meters
        c(CliffHeight<135)=0 ;
        % and set maximum at at cliff height equalt to 450m
        %cMax=fI*450.^(7.2) ;
        %c(c>cMax)=cMax ;

        %     c(c>UserVar.CalvingRateMax)=UserVar.CalvingRateMax ; % set an upper max
        dcddphidx=0;
        dcddphidy=0;

        % Plotting calving law
        %         CH=1:350;
        %         fI=3.2e-17*365.25 ;
        %         c=fI*CF.F.h.^(7.2) ;
        %         c(CH<135)=0;
        %         figure ;
        %         plot(CH,c/1000,LineWidth=2) ;
        %         hold on
        %         plot([135 135],[0 25],'r--')
        %         xlabel("Cliff Height (m)",Interpreter="latex") ;
        %         ylabel("Calving Rate (km/yr)",Interpreter="latex") ;
        %         title("Anna Crawford et al, 2021, T=-20 C $\alpha$=7.2",Interpreter="latex")
        %

        % I=y<-660e3; c(I)=0;
        I=F.y<-550e3; c(I)=0;   % This is what Mathieu uses
        I=F.x<-1700e3; c(I)=0;

        % No calving when B>S
        % Limit calving to 2 km/day
        % cMax=2000*265.25 ; 
        % c(c>cMax)=cMax;   % This makes quite an impact, but also creates a 'spotted' c field 


    case "-NV-"   % Calving is a function of the normal velocity component to the calving front
        % That is:  c = f(v_n)   where c is the calving rate and v_n the ice velocity normal to the calving front
        % v_n = v \cdot n  , where n is a normal to the calving front
        % Using the level set \varphi, we have n=\nabla \varphi / norm(\varphi)

        factor=UserVar.CalvingLaw.Factor ;
        % c=factor*(F.ub.*nx+F.vb.*ny);

        N=sqrt(dphidx.*dphidx+dphidy.*dphidy+eps);

        Vn=-(F.ub.*dphidx+F.vb.*dphidy)./N;  % Note the sign convention 
        c=factor*Vn;

        % c=-factor*(F.ub.*dfdx+F.vb.*dfdy)./N;

        dcddphidx=-factor.* (F.ub./N - dphidx.*(F.ub.*dphidx+F.vb.*dphidy)./( (dphidx.*dphidx+dphidy.*dphidy+eps).^(3/2))  ) ;
        dcddphidy=-factor.* (F.vb./N - dphidy.*(F.ub.*dphidx+F.vb.*dphidy)./( (dphidx.*dphidx+dphidy.*dphidy+eps).^(3/2))  ) ;

    case "-RR-"    % Prescribing retreat rate

       

        factor=1 ; 
        N=sqrt(dphidx.*dphidx+dphidy.*dphidy+eps);
        Vn=-(F.ub.*dphidx+F.vb.*dphidy)./N;  % Normal velocity component

        RetreatRate=Vn*0;
        RetreatRate(Vn>1000)=1000;  % Retreat rate equal to 1000 m/yr where Vn > 1000 m/yr 

        c=Vn+RetreatRate ;

        % c=-factor*(F.ub.*dfdx+F.vb.*dfdy)./N;

        dcddphidx=-factor.* (F.ub./N - dphidx.*(F.ub.*dphidx+F.vb.*dphidy)./( (dphidx.*dphidx+dphidy.*dphidy+eps).^(3/2))  ) ;
        dcddphidy=-factor.* (F.vb./N - dphidy.*(F.ub.*dphidx+F.vb.*dphidy)./( (dphidx.*dphidx+dphidy.*dphidy+eps).^(3/2))  ) ;


    case "-hqk-"

        % this is the inverse-thickness calving law used in the ice-shelf 1D verification experiment
        q=-2;

        k=86320694.4400036;
        c=k*F.h.^q;

        dcddphidx=0;
        dcddphidy=0;

    case "-Fqk-"  % c = k F^q

        % power law as a function of cliff height (freeboard f)

        % The DeConto cliff law seems to have been:
        % c=3000 m/yr for F>=100
        % c=0 for F < 80 
        %
        % A linear interpolation between those two values gives
        %
        %  c(F)= k F
        %
        %  c(100) = k 100 m = 3000 m/yr  -> k= 3000 (m/yr) /100 m  = 30 1/yr 


        CliffHeight=min((F.s-F.S),F.h) ;

        k=UserVar.CalvingLaw.Fqk.k ;
        q=UserVar.CalvingLaw.Fqk.q ;


        c=k.*CliffHeight.^q; 

        c(F<UserVar.CalvingLaw.Fqk.Fmin)= UserVar.CalvingLaw.Fqk.cmin;
        c(F>UserVar.CalvingLaw.Fqk.Fmax)= UserVar.CalvingLaw.Fqk.cmax;

        % Note:  It is presumably best to set
        % k= UserVar.CalvingLaw.Fqk.cmax/(UserVar.CalvingLaw.Fqk.Fmax)^q ;
        % otherwise one could have the situation where c>UserVar.CalvingLaw.Fqk.cmax
        % for F<UserVar.CalvingLaw.Fqk.Fmax


        dcddphidx=0;
        dcddphidy=0;

    case "-DP-"  % Robert DeConto and David Pollard
        
        % Units: meters, year
        rhoice=917;
        CliffHeight=min((F.s-F.S),F.h).*F.rho./rhoice;  % OK, so it is a bit unclear what the "cliff height" should be
                                                      % But presumably the idea is that the CliffHeight is a proxy
                                                      % for the stresses at the calving front, so it appears likley that it should
                                                      % involve rho*g*CliffHeight
                                                      % For this to be comparable with values in the litterature
                                                      % this will most likely need to be adjusted to water equivalent height.


        k1=-12000 ; k2=150 ;
        c=k1+k2*CliffHeight;  % c(CliffHeight=80 m)=0 and c(CliffHeight=100 m) = 3000 m/yr

        c(CliffHeight<80)= 0;
        c(CliffHeight>100)=3000;
        
        dcddphidx=0;
        dcddphidy=0;

    otherwise

        error("CaseNotFound")

end