function  [c,dcddphidx,dcddphidy]=DefineCalvingAtIntegrationPoints(UserVar,CtrlVar,dphidx,dphidy,F) 

%%
%
% Defines calving at integration points. 
%
%  phi is the level set function
%
%  dcddphidx is   dc/d (dphidx) = \frac{ d c]{d (dphidx)}, ie it is the derivative of c with respect to dphidx
%  where dphidx is in turn d(phi)/dx              
%
% Note: F is here provided at the integration points, i.e. this is not the usual F variable provided at nodes!
%       Not all the usual fields of F are available. The fields available include:
%       u,v,h,s,S,rho,exx,exy,eyy.
%
% Also note that F.x,F.y are here the (x,y) coordinates of the integration points, ie not the (x,y) nodal coordinates. 
% 
% The level-set option must be activated by setting
%
%  CtrlVar.LevelSetMethod=1; 
%
% in DefineInitialInputs.m
%
% If the calving rate is a function of the gradients of the level-set, the calving rate must be defined at the element
% integration points using this function as:
%
%
%   [c,dcddphidx,dcddphidy]=DefineCalvingAtIntegrationPoints(UserVar,CtrlVar,dphidx,dphidy,F) 
%
% which provies dphi/dx and dphi/dy where phi is the level-set function
%
% This option is activated by the usuer by setting
%
%
%   CtrlVar.CalvingLaw.Evaluation="-int-" ; 
%
% in DefineInitialInputs.m
%
% The default is to prescribe the calving rate at the nodes, i.e. by default we have
%
%    CtrlVar.CalvingLaw.Evaluation="-node-" ; 
%
% Prescribing the calving rate at the integration points is, for example, required if the calving rate is a function of velocities normal to the calving front.
%
% The user must then also return derivatives of the calving rate, c, with respect to x and y derivatives of the level set,
% ie 
% 
% $$\frac{dc}{d (d \phi/dx)}$$
%
% $$\frac{dc}{d (d \phi/dy)}$$
% 
% More details are provided in the UaCompendium
%
%%

narginchk(5,5)



isPlotStrainRates=false;

if isPlotStrainRates

    FigStrainRates=FindOrCreateFigure("strain rates and velocities") ; clf(FigStrainRates) ;
    QuiverColorGHG(F.x/1000,F.y/1000,F.ub,F.vb) ;
    Scale=2;
    hold on ; PlotTensor(F.x/1000,F.y/1000,F.exx,F.exy,F.eyy,Scale) ; axis equal

end



switch UserVar.CalvingLaw.Type

    case "-AC-"

        %CliffHeight=min((F.s-F.S),F.h) ;  % note: this does not account for density
        
        CliffHeight=min((F.s-F.S),F.h).*F.rho./1000; % guessing this is what the authors intended,
                                                     % but maybe this should be F.rho/920 assuming the authors used 920.

        fI=3.2e-17*365.25 ; c=fI*CliffHeight.^(7.2) ;
        % Now set calving rate to zero for cliff less than 135meters
        c(CliffHeight<135)=0 ;
        % and set maximum at at cliff height equalt to 450m
        cMax=fI*450.^(7.2) ;
        c(c>cMax)=cMax ;

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
        % I=F.y<-550e3; c(I)=0;   % This is what Mathieu uses, but this seems to create problems with LSF
        % I=F.x<-1700e3; c(I)=0;

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

        CliffHeight=min((F.s-F.S),F.h).*F.rho./1000;  % OK, so it is a bit unclear what the "cliff height" should be
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