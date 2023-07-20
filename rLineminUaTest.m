function [gammamin,rmin,du,dv,dh,dl,BackTrackInfo,rForce,rWork,D2] = rLineminUaTest(CtrlVar,UserVar,func,r0,r1,K,L,du0,dv0,dh0,dl0,dJdu,dJdv,dJdh,dJdl,Normalisation,M)

%%

%  [ K  L ]  [dx0]  = [ dJdx ]
%  [ L' 0 ]  [dl0]    [ dJdl ]
%


% func=@(gamma,Du,Dv,Dl) CalcCostFunctionNR(UserVar,RunInfo,CtrlVar,MUA,gamma,F,fext0,L,l,cuv,Du,Dv,Dl) ;
%%

% CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ;

%%

rminNewton=nan     ; rminDescent=nan     ; rminCauchy=nan     ;  rminCauchyM=nan;  rCN=nan;
gammaminNewton=nan ; gammaminDescent=nan ; gammaminCauchy=nan ;  gammaminCauchyM=nan ; gammaCauchyM=nan ; gammaminCN=nan ; 
normNewton=nan ; normCauchy=nan ; normCN=nan; 
%% these will be the returns if no min is found
gammamin=0 ; rmin=r0 ;   du=du0*0 ; dv=dv0*0 ; dh=dh0*0 ; dl=dl0*0 ;
rForce=r0 ; rWork=nan ; D2=nan ;

%%
NoReduction=true; 

BestMethod="" ;

% CtrlVar.rLineMinUa="-Newton-Steepest Descent-Steepest Descent Mass-Cauchy-" ;
% CtrlVar.rLineMinUa="-Newton-Steepest Descent Mass-" ;

CtrlVar.rLineMinUa="-Auto-" ;  % First do the Newton step, ie H n = -g and evaluated r^2(x+n)/r(x)
                               % If reduction not sufficient, do Newton-direction backtracking, ie min r^2(x+ gamma n)
                               % If best backtracking step too small, ie gammaminNewton < gammaMinNewtonAccepted, do the M-Cauchy step, is solve M c = -g and evaluate r^2(x+c)/r(x) 
                               % If full Cauchy step not OK, then do Cauchy backtrackiing 
                               % If full Cauchy step OK, do minimisation along c to n direction (dog-leg)


% CtrlVar.rLineMinUa="-Newton-Steepest Descent-Steepest Descent Mass-Plot Quad Approximations-" ;
% CtrlVar.rLineMinUa="-Newton Step-Cauchy M-step-Cauchy M to Newton-Steepest Descent-Plot Quad Approximations-" ;
% CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ;

%%

if isempty(dh0)
    Variables="-uvl-";
else
    Variables="-uvhl-";
end


if contains(CtrlVar.rLineMinUa,"-Newton Step-")  || contains(CtrlVar.rLineMinUa,"-Auto-")

    % Newton

    % Newton direction is [du0;dv0;dl0] ;

    if Variables=="-uvl-"
        rNewtonFunc=@(gamma) func(gamma,du0,dv0,dl0) ;
        NewtonPoint=[du0;dv0;dl0];
    elseif Variables=="-uvhl-"
        rNewtonFunc=@(gamma) func(gamma,du0,dv0,dh0,dl0) ;
        NewtonPoint=[du0;dv0;dh0;dl0];
    else
        error("variables?")
    end

    if isnan(r0) || isempty(r0)
        r0=rNewtonFunc(0);
        rmin=r0 ;
    end
    if isnan(r1) || isempty(r1)
        r1=rNewtonFunc(1);
    end
    r2Newton=r1 ; % plotting purposes
    NewtonSlope0=-2*r0 ;

    CtrlVar.uvMinimisationQuantity="Force Residuals" ;
    CtrlVar.BacktrackingGammaMin=0.001 ; 
    CtrlVar.BacktracFigName="Line Search in Newton Direction" ;
    [gammaminNewton,rminNewton,BackTrackInfo]=BackTracking(NewtonSlope0,1,r0,r1,rNewtonFunc,CtrlVar);
    NewtonPointUpdated=gammaminNewton*NewtonPoint;
    % Basically whenever the full Newton step is not accepted in the linesearch based on teh Armijo's criterion,
    % calculate the Cauchy-M point

    gammaminNewtonAccepted=0.75 ; rminNewtonRatioAccepted=0.5 ; 
    rminNewtonRatio=rminNewton/r0; 
    if  contains(CtrlVar.rLineMinUa,"-Auto-") &&  gammaminNewton < gammaminNewtonAccepted  && rminNewtonRatio > rminNewtonRatioAccepted

        fprintf(' rLineminUa: Newton step is %f and smaller that %f  \n ',gammaminNewton,gammaminNewtonAccepted) ;
        fprintf(' rLineminUa: Will now try using Cauchy step Steepest Descent were the mass matrix replaces the Hessian \n ') ;
        CtrlVar.rLineMinUa="-Auto-Cauchy M-step-" ;

    end

    normNewton=norm([du0;dv0;dh0;dl0]) ;

    if rminNewton < rmin

        NoReduction=false; 

        du=gammaminNewton*du0;
        dv=gammaminNewton*dv0;
        dh=gammaminNewton*dh0;
        dl=gammaminNewton*dl0;

        gammamin=gammaminNewton;
        rmin=rminNewton;


        [rTest,~,~,rForce,rWork,D2]=rNewtonFunc(gammamin);
        BackTrackInfo.Direction="N " ;
        BestMethod="Newton" ;



        if abs(rTest-rmin)>1000*eps

            error("rLineMinUa:inconsistent","r not having the value expected")

        end

    end

end

if contains(CtrlVar.rLineMinUa,"-Cauchy M-step-") || contains(CtrlVar.rLineMinUa,"-Cauchy M to Newton-") || contains(CtrlVar.rLineMinUa,"-Steepest Descent-")

    if Variables=="-uvl-"
        R=[dJdu;dJdv;dJdl];
        Ruv=[dJdu;dJdv];
        sol0=[du0;dv0];  
    elseif Variables=="-uvhl-"
        R=[dJdu;dJdv;dJdh;dJdl];
        Ruvh=[dJdu;dJdv;dJdh];
        sol0=[du0;dv0;dh0];  
    else
        error(" Variables? ")
    end


    [nL,mL]=size(L);
    L0=sparse(nL,nL);
    H=[K L' ;L L0] ;

end

%% M-Cauchy or the M-Steepest Descent, effectivily as if using the Mass matrix as the Hessian
if contains(CtrlVar.rLineMinUa,"-Cauchy M-step-")
    
    % Du=dx(1:nM) ; Dv=dx(nM+1:2*nM) ; Dl=dx(2*nM+1:end);
    % Replace:
    %
    %  K = [Kuu Kuv ]
    %      [Kuv Kvv [
    %
    % with
    %
    %   K=[M 0 ]
    %     [0 M ]
    %
    %
    [nM,mM]=size(M);
    O=sparse(nM,mM);

    if Variables=="-uvl-"

        Mblock=[M O ; O M ];
        [sol,DlM]=solveKApeSymmetric(Mblock,L,Ruv,dJdl,sol0,dl0,CtrlVar);
        DuM=sol(1:nM) ; DvM=sol(nM+1:2*nM); DhM=[] ; 
        rMassFunc=@(gamma) func(gamma,DuM,DvM,DlM) ;
        sM=[DuM;DvM;DlM];

    elseif Variables=="-uvhl-"


        Mblock=[M O O ; O M O ; O O M];
        [sol,DlM]=solveKApeSymmetric(Mblock,L,Ruvh,dJdl,sol0,dl0,CtrlVar);
        DuM=sol(1:nM) ; DvM=sol(nM+1:2*nM);  DhM=sol(2*nM+1:3*nM); 
        rMassFunc=@(gamma) func(gamma,DuM,DvM,DhM,DlM) ;
        sM=[DuM;DvM;DhM;DlM];


    end

    
   % To do:  define a temp variable: RH=R'*H ; 
    CauchyMSlope0=(-2*R'*H*sM)/Normalisation;
    % CauchyMSlope0=-(sM'*(H'*R)+(R'*H)*sM)/Normalisation;

    gammaCauchyM=(R'*H*sM)/((H*sM)'*(H*sM)) ;
    CauchyMPoint=gammaCauchyM*sM ;
    normCauchy=norm(CauchyMPoint) ;
    rD0=r0 ;


    if CauchyMSlope0 > 0
        CauchyMSlope0=-CauchyMSlope0;
        gammaCauchyM = -0.1 *rD0/CauchyMSlope0 ;  % initial step size
    end

    rCauchyM=rMassFunc(gammaCauchyM);
    CtrlVar.uvMinimisationQuantity="Force Residuals" ;  
    CtrlVar.BacktracFigName="Line Search in Cauchy M direction" ;

    CtrlVar.LineSearchAllowedToUseExtrapolation=true;
    CtrlVar.NewtonAcceptRatio=0.999; % use the Armijo's condition based on slope0
                    
    [gammaminCauchyM,rminCauchyM,BackTrackInfoCauchyM]=BackTracking(CauchyMSlope0,gammaCauchyM,r0,rCauchyM,rMassFunc,CtrlVar);
    CauchyMPointUpdated=gammaminCauchyM*sM ;
    r2MD=rminCauchyM ; % plotting purposes
  


    if rminCauchyM < rmin     % if this rmin is the smallest so far, update solution

        NoReduction=false; 

        du=gammaminCauchyM*DuM;
        dv=gammaminCauchyM*DvM;
        dh=gammaminCauchyM*DhM;
        dl=gammaminCauchyM*DlM;
        gammamin=gammaminCauchyM;
        rmin=rminCauchyM;
        BackTrackInfo=BackTrackInfoCauchyM;
        BackTrackInfo.Direction="MD" ;
        BestMethod="Mass Steepest Descent" ;
        [rTest,~,~,rForce,rWork,D2]=rMassFunc(gammamin);
        gammamin=gammaminCauchyM/gammaCauchyM ;  % on return, normalize this with the min of the quad model
        
    end

    if rminCauchyM < r0

        NoReduction=false; 

        fprintf(' rLineminUa: M-Cauchy point led to a reduction with respect to r0, with rminCauchyM/r0=%f \n ',rminCauchyM/r0) ;

        if  contains(CtrlVar.rLineMinUa,"-Auto-")

            fprintf(' rLineminUa: Will now do the rest of the dogleg and go from Cauchy to Newton \n ') ;

            CtrlVar.rLineMinUa="-Auto-Cauchy M to Newton-" ;
            % If Newton step is not accepted, try Steepest Descent.

        end

        if  contains(CtrlVar.rLineMinUa,"-Auto-")

            if gammaminCauchyM > gammaCauchyM / 2   % if minimum was found close to the Cauchy point, try Cauch M to Newton direction as well

                CtrlVar.rLineMinUa="-Auto-Cauchy M to Newton-" ;

            % else


                dx=sM ;                % setting direction the be the that of M-modified steepest Descent
                Q=@(g) (R-g*(H*dx))'*(R-g*(H*dx))/Normalisation ;


                rminQCauchyM=Q(gammaCauchyM);
                EstimatedReduction=r0-rminQCauchyM;
                ActualReductoni=r0-rminCauchyM;
                ratioEstimated2Actual=EstimatedReduction/ActualReductoni;
                fprintf("ratio between estimated and realized reduction in MCauchy step is %f \n",ratioEstimated2Actual)

                if ratioEstimated2Actual > 0.8 && ratioEstimated2Actual < 1.8
                    % Q model is good, so try going from here towards the Newton point
                    CtrlVar.rLineMinUa="-Auto-Cauchy M to Newton-" ;
                end

            end
        end

    else

        CtrlVar.rLineMinUa="-Auto-Steepest Descent-";

    end
end

% if  contains(CtrlVar.rLineMinUa,"-Auto-")  && rmin/r0 > 0.9
% 
%     % If reduction in the Newton and the M-Cauchy step was too small,
%     % try steepest Descent as well
%     CtrlVar.rLineMinUa=CtrlVar.rLineMinUa+"-Steepest Descent-";
% 
% end

%% I-Cauchy, or the Steepest Descent, effectivily as if using the Unit matrix as the Hessian
if contains(CtrlVar.rLineMinUa,"-Steepest Descent-")
    [nM,mM]=size(M);
    % Must still solve a system because I need to preserve the BCs
    if Variables=="-uvl-"

        I2n=speye(2*nM,2*nM) ;
        Ruv=[dJdu;dJdv];
        sol0=[du0;dv0];  
        [sol,dlD]=solveKApeSymmetric(I2n,L,Ruv,dJdl,sol0,dl0,CtrlVar);
        duD=sol(1:nM) ; dvD=sol(nM+1:2*nM); dhD=[] ;
        rDescentFunc=@(gamma) func(gamma,duD,dvD,dlD) ;
        sD=[duD;dvD;dlD];
    

    elseif Variables=="-uvhl-"

        I3n=speye(3*nM,3*nM) ;
      
        Ruvh=[dJdu;dJdv;dJdh];
        sol0=[du0;dv0;dh0];  
        [sol,dlD]=solveKApeSymmetric(I3n,L,Ruvh,dJdl,sol0,dl0,CtrlVar);
        duD=sol(1:nM) ; dvD=sol(nM+1:2*nM);  dhD=sol(2*nM+1:3*nM);
        rDescentFunc=@(gamma) func(gamma,duD,dvD,dhD,dlD) ;
        sD=[duD;dvD;dhD;dlD];
        

    end


    slope0Descent=-2*R'*H*sD/Normalisation ;

    if slope0Descent > 0
        slope0Descent = - slope0Descent;
    end

    gammaminCauchy=0.5*sD'*(H'+H)*sD/(sD'*(H'*H)*sD) ;
    gammaCauchyD=(R'*H*sD)/((H*sD)'*(H*sD)) ;
    
    if gammaCauchyD > 0
        b=gammaCauchyD ;
    else
        b = -0.1 *r0/slope0Descent ;  % initial step size
    end


    gamma=b ; rb=rDescentFunc(gamma);

    % CtrlVar.InfoLevelBackTrack=1000;
    CtrlVar.BacktracFigName="Steepest Descent" ;
    CtrlVar.BacktrackingGammaMin=gammaminCauchy/1000; 
    CtrlVar.LineSearchAllowedToUseExtrapolation=true;
    CtrlVar.NewtonAcceptRatio=0.999; % primarirly use the Armijo's condition based on slope0

    [gammaminDescent,rminDescent,BackTrackInfoSteepest]=BackTracking(slope0Descent,b,r0,rb,rDescentFunc,CtrlVar);

    if rminDescent < rminCauchy
        CauchyMPointUpdated=gammaminDescent*sD;
        r2MD=rminDescent ; % plotting purposes
    end

 

    if rminDescent < rmin
        NoReduction=false;
        du=gammaminDescent*duD;
        dv=gammaminDescent*dvD;
        dh=gammaminDescent*dhD;
        dl=gammaminDescent*dlD;
        gammamin=gammaminDescent;
        rmin=rminDescent;
        BackTrackInfo=BackTrackInfoSteepest;
        BackTrackInfo.Direction="SD" ;
        BestMethod="Steepest Descent" ;
        [rTest,~,~,rForce]=rDescentFunc(gammamin);
        rWork=nan ; D2=nan ;  % those have no meaning for the steepest Descent direction
        gammamin=gammamin/gammaminCauchy ; % on return, normalize this with the min of the quad model
    end

end



%% HR-Cauchy, Steepest as 2 H' R

if contains(CtrlVar.rLineMinUa,"-Steepest Descent HR-")
    [nM,mM]=size(M);
   
    if Variables=="-uvl-"

        I2n=speye(2*nM,2*nM) ;
        Ruv=[dJdu;dJdv];
        HR=-2*K'*Ruv ;
        sol0=[du0;dv0];
        [sol,dlD]=solveKApeSymmetric(I2n,L,HR,dJdl,sol0,dl0,CtrlVar);
        duD=sol(1:nM) ; dvD=sol(nM+1:2*nM); dhD=[] ;
        rHRFunc=@(gamma) func(gamma,duD,dvD,dlD) ;
        sHR=[duD;dvD;dlD];

    elseif Variables=="-uvhl-"


        I3n=speye(3*nM,3*nM) ;
        Ruvh=[dJdu;dJdv;dJdh];
        H=[K L' ;L L0] ;

        HR=-2*K'*Ruvh ;

        if ~isempty(L)
            frhs=-Ruvh-L'*luvh;
            grhs=cuvh-L*[F1.ub;F1.vb;F1.h];
        else
            frhs=-Ruvh;
            grhs=[];
        end

        sol0=[du0;dv0;dh0];
        [sol,dlHR]=solveKApeSymmetric(I3n,L,HR,dJdl,sol0,dl0,CtrlVar);

        duHR=sol(1:nM) ; dvHR=sol(nM+1:2*nM);  dhHR=sol(2*nM+1:3*nM);
        rHRFunc=@(gamma) func(gamma,duHR,dvHR,dhHR,dlHR) ;
        sHR=[duHR;dvHR;dhHR;dlHR];

    end

    R=[dJdu;dJdv;dJdh;dJdl];
    sHR=2*H'*R ; 
    duHR=sHR(1:nM) ; dvHR=sHR(nM+1:2*nM);  dhHR=sHR(2*nM+1:3*nM); dlHR=sHR(3*nM+1:end) ;
    rHRFunc=@(gamma) func(gamma,duHR,dvHR,dhHR,dlHR) ;

    slope0HR=-2*R'*H*sHR/Normalisation ;
    if slope0HR > 0
        slope0HR = - slope0HR;
    end

    gammaCauchyHR=R'*H*sHR/(sHR'*(H'*H)*sHR) ;

    
    if gammaCauchyHR > 0
        b=gammaCauchyHR ;
    else
        b = -0.1 *r0/slope0HR ;  % initial step size
    end


    gamma=b ; rb=rHRFunc(gamma);

    % CtrlVar.InfoLevelBackTrack=1000;
    CtrlVar.BacktracFigName="Steepest Descent HR" ;
    CtrlVar.BacktrackingGammaMin=gammaCauchyHR/1000; 
    CtrlVar.LineSearchAllowedToUseExtrapolation=true;
    CtrlVar.NewtonAcceptRatio=0.999; % primarirly use the Armijo's condition based on slope0

    [gammaminCauchyHR,rminHR,BackTrackInfoHR]=BackTracking(slope0Descent,b,r0,rb,rHRFunc,CtrlVar);

    if rminHR < rminCauchyM  % This is for the C to N vector
        CauchyMPointUpdated=gammaminCauchyHR*sHR;
        r2HR=rminHR ; % plotting purposes
    end

 

    if rminHR < rmin
        NoReduction=false;
        du=gammaminCauchyHR*duHR;
        dv=gammaminCauchyHR*dvHR;
        dh=gammaminCauchyHR*dhHR;
        dl=gammaminCauchyHR*dlHR;
        gammamin=gammaminCauchyHR;
        rmin=rminHR;
        BackTrackInfo=BackTrackInfoSteepest;
        BackTrackInfo.Direction="HR" ;
        BestMethod="Steepest Descent HR" ;
        [rTest,~,~,rForce]=rHRFunc(gammamin);
        rWork=nan ; D2=nan ;  % those have no meaning for the steepest Descent direction
        gammamin=gammaminCauchyHR/gammaCauchyHR ; % on return, normalize this with the min of the quad model
    end

end




if gammaminNewton < 0.01  
    CtrlVar.rLineMinUa=replace(CtrlVar.rLineMinUa,"Cauchy M to Newton","") ;   
    CtrlVar.rLineMinUa=replace(CtrlVar.rLineMinUa,"--","-") ;
end

%% Dogleg, going from the M-Cauchy point towards the Newton point
if contains(CtrlVar.rLineMinUa,"-Cauchy M to Newton-")


    [nM,mM]=size(M);
    TolX=0.02;
    options = optimset('Display','iter','TolX',TolX,'OutputFcn',@outfunFminbnd);
    FunCN=@(step) Cauchy2Newton(func,step,"-uvhl-",nM,CauchyMPointUpdated,NewtonPoint) ;  % remember also to change below in any other calls to Cauchy2Newton! 
    
    [gammaminCN,rCN,exitflag,output]=fminbnd(FunCN,0,1,options) ;
    [stop,Outs]=outfunFminbnd();

    xVector=Outs.x ;  fVector=Outs.f ;
    I=~isnan(xVector)  ; xVector=xVector(I) ; fVector=fVector(I) ;

    if rCN < rmin

        NoReduction=false;   
   
        [r2Test,rForce,rWork,DuC2N,DvC2N,DhC2N,DlC2N]=Cauchy2Newton(func,gammaminCN,"-uvhl-",nM,CauchyMPointUpdated,NewtonPoint) ;
        
        du=DuC2N;
        dv=DvC2N;
        dh=DhC2N;
        dl=DlC2N;
        gammamin=gammaminCN;
        rmin=r2Test;  % should be same as rCN
        BackTrackInfo.iarm=output.iterations ;
        BackTrackInfo.Converged=true; 
        BackTrackInfo.nExtrapolationSteps=0;
        BackTrackInfo.Infovector=nan ; 
        BackTrackInfo.Direction="CN" ;
        BestMethod="Cauchy2Newton" ;

        normCN=norm([DuC2N;DvC2N;DhC2N;DlC2N]) ;

    end


    if CtrlVar.InfoLevelBackTrack >= 1000
        
        
        % r2Newton=rNewtonFunc(1) ;
        % r2MD=rMassFunc(gammaminCauchyM) ;


        figCN=FindOrCreateFigure(" Cauchy to Newton ") ; clf(figCN) ;
        xVector=xVector(:) ; fVector=fVector(:) ; 
        xVector=[xVector;0;1]; fVector=[fVector;r2MD;r2Newton];
        [~,I]=sort(xVector)  ; xVector=xVector(I) ; fVector=fVector(I) ;

        plot(xVector,fVector,'o-r') ;
        hold on ;
        plot(0,r2MD,"*m")
        plot(1,r2Newton,"*b")
        text(0,r2MD,"   C",HorizontalAlignment="left")
        text(1,r2Newton,"N   ",HorizontalAlignment="right")
        yline(r0,"--k","$r_0$",interpreter="latex")
        yline(rminNewton,"-.k","$r_N$",interpreter="latex")
        plot(gammaminCN,rCN,'o',MarkerFaceColor="b",MarkerSize=10)
        xlabel("Distance along the Cauchy-to-Newton path",interpreter="latex")
        ylabel("$r^2$",interpreter="latex")
        legend("$r^2$","$r^2_C$","$r^2_N$","$r^2_0$","$\min r^2_N$","$\min r^2_{CM}$",interpreter="latex",location="best")
        title("From the Cauchy-Point to the Newton-Point")
        subtitle(sprintf("t=%f   dt=%f",CtrlVar.time,CtrlVar.dt),Interpreter="latex")
        % fig=gcf ; exportgraphics(fig,"ExampleCaucyToNewtonPath.pdf")
        drawnow
    end



end


if contains(CtrlVar.rLineMinUa,"-Plot Quad Approximations-")

    nPoints=13 ; 
    if Variables=="-uvl-"
        rDescentFunc=@(gamma) func(gamma,dJdu,dJdv,dJdl) ;
    elseif Variables=="-uvhl-"
        rDescentFunc=@(gamma) func(gamma,dJdu,dJdv,dJdh,dJdl) ;
    else
        error("variables?")
    end

    if Variables=="-uvl-"
        % Quad approximation in Newton direction
        dx=[du0;dv0;dl0] ;      % Newton direction
        sM=[DuM;DvM;DlM];       % M direction
        sD=[duD;dvD;dlD];       % D direction
        R=[dJdu;dJdv;dJdl];     % Residual vector
    elseif Variables=="-uvhl-"
        dx=[du0;dv0;dh0;dl0] ;     % Newton direction
        sM=[DuM;DvM;DhM;DlM];      % M direction
        sD=[duD;dvD;dhD;dlD];      % D direction
        R=[dJdu;dJdv;dJdl];        % Residual vector
    end

    % Newton Q model
    Quad=@(g) (R-g*(H*dx))'*(R-g*(H*dx))/Normalisation ;

    Q1=Quad(1) ; 
    gVector=linspace(0,1.2,nPoints)' ;
    QVector=gVector*0+nan;
    rVector=gVector*0+nan;

    for I=1:numel(gVector)

        QVector(I)=Quad(gVector(I));
        rVector(I)=rNewtonFunc(gVector(I)) ;
    end


    fig=FindOrCreateFigure("Q and r (Newton)") ; clf(fig) ;

    hold on
    plot(gVector,rVector,"o-b")
    plot(gVector,QVector,"*-r")
    plot(gammaminNewton,rminNewton,'o',MarkerFaceColor="b",MarkerSize=10)

    legend("r^2","Quad approximation in Newton direction")

    % determine if Q good at full step
    rQN=Quad(1)/rNewtonFunc(1) ;

    fprintf(" Newton  Quad model \n ")
    [gVector QVector rVector QVector./rVector]
    if gammaminNewton > 0.5

        fprintf("Newton step size %f > 0.5. \n",gammaminNewton)
        NewtonModelGood=true;
    else
        NewtonModelGood=false;
    end

    %% Quad approximation in the direction of steepest Descent
    % setting direction the be the that of steepest Descent
    Q=@(g) (R-g*(H*sD))'*(R-g*(H*sD))/Normalisation ;

   % gammaminCauchy=0.5*sD'*(H+H')*sD/(sD'*(H'*H)*sD) ;

    gVector=linspace(0,2*gammaminCauchy,nPoints)' ;
    gVector(2)=gammaminCauchy/1000;

    QVector=gVector*0+nan;
    rVector=gVector*0+nan;
    for I=1:numel(gVector)

        QVector(I)=Q(gVector(I));
        rVector(I)=rDescentFunc(gVector(I)) ;
    end

    fig=FindOrCreateFigure("Q and r steepest Descent") ; clf(fig) ;

    hold on
    plot(gVector,rVector,"o-b")
    plot(gVector,QVector,"*-r")
    legend("r^2","Quad approximation in steepest Descent direction")
    plot(gammaminDescent,rminDescent,'o',MarkerFaceColor="b",MarkerSize=10)

  %% Quad approximation in the direction of M-modified steepest Descent
    % setting direction the be the that of M-modified steepest Descent
    Q=@(g) (R-g*(H*sM))'*(R-g*(H*sM))/Normalisation ;

    

    gVector=linspace(0,2*gammaminCauchyM,nPoints)' ;
    gVector(2)=gammaminCauchyM/1000;

    QVector=gVector*0+nan;
    rVector=gVector*0+nan;
    for I=1:numel(gVector)

        QVector(I)=Q(gVector(I));
        rVector(I)=rMassFunc(gVector(I)) ;

    end

    rQM=Q(1)/rMassFunc(1) ;

    if rQM> 0.8 && rQN < 1.2

        fprintf("Mass Q model within 20%% with rQN=Q(1)/N(1)=%f \n",rQN)
        MassModelGood=true;
    else
        MassModelGood=false;
    end
    fprintf(" CauchyM  Quad model \n ")
    [gVector QVector rVector QVector./rVector]
    fig=FindOrCreateFigure("Q and r M-modified steepest Descent") ; clf(fig) ;

    hold on
    plot(gVector,rVector,"o-b")
    plot(gVector,QVector,"*-r")
    
    plot(gammaminCauchyM,rminCauchyM,'o',MarkerFaceColor="b",MarkerSize=10)
    plot(gammaCauchyM,rCauchyM,'o',MarkerFaceColor="r",MarkerSize=5)
    legend("r^2","Quad approximation in M-modified steepest Descent direction","Estimated min","Initial Guess (MCauchy Point)")

    %% Plot r2 along Cauchy to Newton direction


    

    gVector=linspace(0,1,nPoints)' ;
    [nM,mM]=size(M);
    r2=nan(nPoints,1) ;

    for I=1:nPoints
        r2(I)=Cauchy2Newton(func,gVector(I),"-uvhl-",nM,CauchyMPointUpdated,NewtonPoint) ;
    end

    options = optimset('Display','iter','TolX',0.02,'OutputFcn',@outfunFminbnd);
    FunCN=@(step) Cauchy2Newton(func,step,"-uvhl-",nM,CauchyMPoint,NewtonPoint) ;
    [gammaminCN,rCN,exitflag,output]=fminbnd(FunCN,0,1,options) ;
    [stop,Outs]=outfunFminbnd();
    xVector=Outs.x ;  fVector=Outs.f ;
    I=~isnan(xVector)  ; xVector=xVector(I) ; fVector=fVector(I) ;

    if rCN < rmin

        [r2Test,rForce,rWork,DuC2N,DvC2N,DhC2N,DlC2N]=Cauchy2Newton(func,gammaminCN,"-uvhl-",nM,CauchyMPointUpdated,NewtonPoint) ;
        
        du=DuC2N;
        dv=DvC2N;
        dh=DhC2N;
        dl=DlC2N;
        gammamin=gammaminCN;
        rmin=r2Test;  % should be same as rCN
        BackTrackInfo.iarm=output.iterations ;
        BackTrackInfo.Converged=true; 
        BackTrackInfo.nExtrapolationSteps=0;
        BackTrackInfo.Infovector=nan ; 
        BackTrackInfo.Direction="CN" ;
        BestMethod="Cauchy2Newton" ;
    end

    r2Newton=rNewtonFunc(1) ;
    r2MD=rMassFunc(gammaminCauchyM) ;
    figCN=FindOrCreateFigure(" Cauchy to Newton ") ; clf(figCN) ;
    
    
    xVector=[xVector(:);gVector(:)] ; fVector=[fVector(:);r2(:)] ; 
    [~,I]=sort(xVector)  ; xVector=xVector(I) ; fVector=fVector(I) ; 

    plot(xVector,fVector,'o-r') ; 
    hold on ; 
    plot(0,r2MD,"*m")
    plot(1,r2Newton,"*b")
    text(0,r2MD,"   C(M)",HorizontalAlignment="left")
    text(1,r2Newton,"N   ",HorizontalAlignment="right")
    yline(r0,"--k","$r_0$",interpreter="latex")
    yline(rminNewton,"-.k","$r_N$",interpreter="latex")
    plot(gammaminCN,rCN,'o',MarkerFaceColor="b",MarkerSize=10)
    xlabel("Distance along the Cauchy-to-Newton path",interpreter="latex")
    ylabel("$r^2$",interpreter="latex")
    legend("$r^2$","$r^2_C$","$r^2_N$","$r^2_0$","$\min r^2_N$","$\min r^2_{CM}$",interpreter="latex",location="best")
    % fig=gcf ; exportgraphics(fig,"ExampleCaucyToNewtonPath.pdf")
    drawnow

    




end

if NoReduction
    BackTrackInfo.Converged = false ;
else
    BackTrackInfo.Converged = true ;
end

%% Summary
if CtrlVar.InfoLevelBackTrack >= 2
    fprintf(" [---------- rLineminUa: \n")
    fprintf("\t r0=%-13.7g \t r1/r0=%-13.7g \t rNewton/r0=%-13.7g \t rminCauchyM/r0=%-13.7g \t rDescent/r0=%-13.7g \t rCN/r0=%-13.7g \n",r0,r1/r0,rminNewton/r0,rminCauchyM/r0,rminDescent/r0,rCN/r0)
    fprintf("\t g0=%-13.8g \t    g1=%-13.7g \t    gNewton=%-13.7g \t             gM=%-13.7g \t    gDescent=%-13.7g \t    gCM=%-13.7g \n",0,1,gammaminNewton,gammaminCauchyM/gammaCauchyM,gammaminDescent,gammaminCN)
    fprintf("\t normNewton=%-13.7g \t normCauchy/normNewton=%-13.7g \t normCN/normNewton=%-13.7g \n ",normNewton,normCauchy/normNewton,normCN/normNewton)
    if BestMethod=="Cauchy2Newton"
        fprintf("\t minC2N/minN=%f   \n",rCN/rminNewton)
    end
    fprintf("\t =================>  Best method is %s  with rmin/r0=%f    <=====================\n",BestMethod,rmin/r0)
    fprintf(" -------------------------] \n")
end
%%


return

end
 

 % 
 % [---------- rLineminUa: 
	%  r0=3.041337e-06  	 r1/r0=1358.289      	 rNewton/r0=0.9940949     	 rM/r0=0.7306787     	 rDescent/r0=0.8374294     	 rCauchy/r0=NaN           	 rCM/r0=0.7052111     
	%  g0=0             	    g1=1             	    gNewton=0.02910391    	    gM=1             	    gDescent=1.90946e-11   	    gCauchy=1.90946e-11   
	%  =================>  Best method is Cauchy2Newton     <=====================
 % 


















