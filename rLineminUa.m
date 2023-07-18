function [gammamin,rmin,du,dv,dh,dl,BackTrackInfo,rForce,rWork,D2] = rLineminUa(CtrlVar,UserVar,func,r0,r1,K,L,du0,dv0,dh0,dl0,dJdu,dJdv,dJdh,dJdl,Normalisation,M)

%%
%
%   r=func(gamma,Du,Dv,Dh,Dl)
%
%  [ K  L ]  [dx0]  = [ dJdx ]
%  [ L' 0 ]  [dl0]    [ dJdl ]
%


% func=@(gamma,Du,Dv,Dl) CalcCostFunctionNR(UserVar,RunInfo,CtrlVar,MUA,gamma,F,fext0,L,l,cuv,Du,Dv,Dl) ;
%%

% CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ;

%%

rminNewton=nan     ; rminCauchyD=nan     ; rminCauchy=nan     ;  rminCauchyM=nan;  rCN=nan;
gammaminNewton=nan ; gammaminCauchyD=nan ; gammaCauchyD=nan ;  gammaminCauchyM=nan ; gammaCauchyM=nan ; gammaminCN=nan ; 
normNewtonStep=nan ; normCauchyM=nan ; normCN=nan; normCauchyD=nan ; normNewton=nan ; 
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


%%

if isempty(dh0)
    Variables="-uvl-";
else
    Variables="-uvhl-";
end

rRatioReductionAccepted=0.65 ; 


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

    gammaminNewtonAccepted=0.75 ; rminNewtonRatioAccepted=rRatioReductionAccepted ;  
    rminNewtonRatio=rminNewton/r0; 
    if  contains(CtrlVar.rLineMinUa,"-Auto-") &&  gammaminNewton < gammaminNewtonAccepted  && rminNewtonRatio > rminNewtonRatioAccepted

        if CtrlVar.InfoLevelNonLinIt >= 10 
            fprintf(' rLineminUa: Newton step is %f and smaller that %f  \n ',gammaminNewton,gammaminNewtonAccepted) ;
            fprintf(' rLineminUa: Will now try using Cauchy step Steepest Descent were the mass matrix replaces the Hessian \n ') ;
        end
        CtrlVar.rLineMinUa="-Auto-Cauchy M-step-" ;

    end

    normNewtonStep=norm([du0;dv0;dh0;dl0]) ;

    if rminNewton < rmin

        NoReduction=false; 

        du=gammaminNewton*du0;
        dv=gammaminNewton*dv0;
        dh=gammaminNewton*dh0;
        dl=gammaminNewton*dl0;

        gammamin=gammaminNewton;
        rmin=rminNewton;
        normNewton=norm([du;dv;dh;dl]) ;

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


    [nL,~]=size(L);
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
        funcCauchyM=@(gamma) func(gamma,DuM,DvM,DlM) ;
        sM=[DuM;DvM;DlM];

    elseif Variables=="-uvhl-"


        Mblock=[M O O ; O M O ; O O M];
        [sol,DlM]=solveKApeSymmetric(Mblock,L,Ruvh,dJdl,sol0,dl0,CtrlVar);
        DuM=sol(1:nM) ; DvM=sol(nM+1:2*nM);  DhM=sol(2*nM+1:3*nM); 
        funcCauchyM=@(gamma) func(gamma,DuM,DvM,DhM,DlM) ;
        sM=[DuM;DvM;DhM;DlM];


    end


    % To do:  define a temp variable: RH=R'*H ;
    CauchyMSlope0=(-2*R'*H*sM)/Normalisation;
    gammaCauchyM=(R'*H*sM)/((H*sM)'*(H*sM)) ;
    CauchyMPoint=gammaCauchyM*sM ;




    if CauchyMSlope0 > 0

        %% OK, this is a breakdown as the direction does not lead to a reduction
        %  Try simply to reverse direction...
        CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ; CtrlVar.doplots=1;
        if Variables=="-uvl-"
            funcCauchyM=@(gamma) func(gamma,-DuM,-DvM,-DlM) ;
            sM=-[DuM;DvM;DlM];
        else
            funcCauchyM=@(gamma) func(gamma,-DuM,-DvM,-DhM,-DlM) ;
            sM=-[DuM;DvM;DhM;DlM];

        end
        CauchyMSlope0=(-2*R'*H*sM)/Normalisation;
        gammaCauchyM=(R'*H*sM)/((H*sM)'*(H*sM)) ;
        CauchyMPoint=gammaCauchyM*sM ;
    end


    if CauchyMSlope0 > 0

        CauchyMSlope0=-CauchyMSlope0;
        gammaCauchyM = -0.1 *r0/CauchyMSlope0 ;  % initial step size
    end

    if gammaCauchyM < 0
        gammaCauchyM = -0.1 *r0/CauchyMSlope0 ;  % initial step size
    end


    rCauchyM=funcCauchyM(gammaCauchyM);
    CtrlVar.uvMinimisationQuantity="Force Residuals" ;  
    CtrlVar.BacktracFigName="Line Search in Cauchy M direction" ;

    CtrlVar.LineSearchAllowedToUseExtrapolation=true;
    CtrlVar.NewtonAcceptRatio=0.9; % use the Armijo's condition based on slope0
    CtrlVar.BacktrackingGammaMin=1e-10 ; 
    [gammaminCauchyM,rminCauchyM,BackTrackInfoCauchyM]=BackTracking(CauchyMSlope0,gammaCauchyM,r0,rCauchyM,funcCauchyM,CtrlVar);
    CauchyPointUpdated=gammaminCauchyM*sM ;
    r2MD=rminCauchyM ; % plotting purposes
  

     normCauchyM=gammaminCauchyM*norm([du;dv;dh;dl]) ;

    if rminCauchyM < rmin     % if this rmin is the smallest so far, update solution

        NoReduction=false; 

        du=gammaminCauchyM*DuM;
        dv=gammaminCauchyM*DvM;
        dh=gammaminCauchyM*DhM;
        dl=gammaminCauchyM*DlM;
        
        rmin=rminCauchyM;
        BackTrackInfo=BackTrackInfoCauchyM;
        BackTrackInfo.Direction="MD" ;
        BestMethod="Mass Steepest Descent" ;
        % gammamin=gammaminCauchyM;
        % [rTest,~,~,rForce,rWork,D2]=funcCauchyM(gammamin);
        rForce=rmin ; rWork=nan ; D2=nan ;
        gammamin=gammaminCauchyM/gammaCauchyM ;  % on return, normalize this with the min of the quad model


    end

    if rminCauchyM < r0

        NoReduction=false;
        if CtrlVar.InfoLevelNonLinIt >= 10
            fprintf(' rLineminUa: M-Cauchy point led to a reduction with respect to r0, with rminCauchyM/r0=%f \n ',rminCauchyM/r0) ;
        end
        if  contains(CtrlVar.rLineMinUa,"-Auto-")

            if CtrlVar.InfoLevelNonLinIt >= 10
                fprintf(' rLineminUa: Will now do the rest of the dogleg and go from Cauchy to Newton \n ') ;
            end
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

                if CtrlVar.InfoLevelNonLinIt >= 10
                    fprintf("ratio between estimated and realized reduction in MCauchy step is %f \n",ratioEstimated2Actual)
                end
                if ratioEstimated2Actual > 0.8 && ratioEstimated2Actual < 1.8
                    % Q model is good, so try going from here towards the Newton point
                    CtrlVar.rLineMinUa="-Auto-Cauchy M to Newton-" ;
                end

            end
        end

    else

        CtrlVar.rLineMinUa="-Auto-Steepest Descent-Cauchy M to Newton-"; 

    end
end


if rmin/r0 < rRatioReductionAccepted

    CtrlVar.rLineMinUa="" ;

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
    [nM,~]=size(M);
    % Must still solve a system because I need to preserve the BCs
    if Variables=="-uvl-"

        I2n=speye(2*nM,2*nM) ;
        Ruv=[dJdu;dJdv];
        sol0=[du0;dv0];  
        [sol,dlD]=solveKApeSymmetric(I2n,L,Ruv,dJdl,sol0,dl0,CtrlVar);
        duD=sol(1:nM) ; dvD=sol(nM+1:2*nM); dhD=[] ;
        funcCauchyD=@(gamma) func(gamma,duD,dvD,dlD) ;
        sD=[duD;dvD;dlD];

    elseif Variables=="-uvhl-"

        I3n=speye(3*nM,3*nM) ;
      
        Ruvh=[dJdu;dJdv;dJdh];
        sol0=[du0;dv0;dh0];  
        [sol,dlD]=solveKApeSymmetric(I3n,L,Ruvh,dJdl,sol0,dl0,CtrlVar);
        duD=sol(1:nM) ; dvD=sol(nM+1:2*nM);  dhD=sol(2*nM+1:3*nM);
        funcCauchyD=@(gamma) func(gamma,duD,dvD,dhD,dlD) ;
        sD=[duD;dvD;dhD;dlD];

    end

    CauchyDSlope0=-2*R'*H*sD/Normalisation ;
    gammaCauchyD=(R'*H*sD)/((H*sD)'*(H*sD)) ;
    % CauchyDPoint=gammaCauchyD*sD;
    

    if CauchyDSlope0 > 0
        CauchyDSlope0 = - CauchyDSlope0;
    end


    if gammaCauchyD < 0
        gammaCauchyD = -0.1 *r0/CauchyDSlope0 ;  % initial step size
    end


    rCauchyD=funcCauchyD(gammaCauchyD);

    % CtrlVar.InfoLevelBackTrack=1000;
    CtrlVar.BacktracFigName="Line Search in Cauchy D direction" ;
    CtrlVar.BacktrackingGammaMin=gammaCauchyD/1000; 
    CtrlVar.LineSearchAllowedToUseExtrapolation=true;
    CtrlVar.NewtonAcceptRatio=0.999; % primarirly use the Armijo's condition based on slope0

    [gammaminCauchyD,rminCauchyD,BackTrackInfoCauchyD]=BackTracking(CauchyDSlope0,gammaCauchyD,r0,rCauchyD,funcCauchyD,CtrlVar);

    if rminCauchyD < rminCauchyM
        CauchyPointUpdated=gammaminCauchyD*sD;
        r2MD=rminCauchyD ; % plotting purposes
    end

 
    normCauchyD=gammaminCauchyD*norm([du;dv;dh;dl]) ;

    if rminCauchyD < rmin
        NoReduction=false;
        du=gammaminCauchyD*duD;
        dv=gammaminCauchyD*dvD;
        dh=gammaminCauchyD*dhD;
        dl=gammaminCauchyD*dlD;
        gammamin=gammaminCauchyD;
        rmin=rminCauchyD;
        BackTrackInfo=BackTrackInfoCauchyD;
        BackTrackInfo.Direction="SD" ;
        BestMethod="Steepest Descent" ;
        [~,~,~,rForce]=funcCauchyD(gammamin);
        rWork=nan ; D2=nan ;  % those have no meaning for the steepest Descent direction
        gammamin=gammamin/gammaCauchyD ; % on return, normalize this with the min of the quad model
     
    end

end


if rmin/r0 < rRatioReductionAccepted

    CtrlVar.rLineMinUa="" ;

end

% 
% if gammaminNewton < 0.01  
%     CtrlVar.rLineMinUa=replace(CtrlVar.rLineMinUa,"Cauchy M to Newton","") ;   
%     CtrlVar.rLineMinUa=replace(CtrlVar.rLineMinUa,"--","-") ;
% end

%% Dogleg, going from the M-Cauchy point towards the Newton point
if contains(CtrlVar.rLineMinUa,"-Cauchy M to Newton-")

    % CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ;
    [nM,mM]=size(M);
    TolX=0.02;


    if CtrlVar.InfoLevelNonLinIt >= 10 
        options = optimset('Display','iter','TolX',TolX,'OutputFcn',@outfunFminbnd);
    else
        options = optimset('Display','off','TolX',TolX,'OutputFcn',@outfunFminbnd);
    end


    FunCN=@(step) Cauchy2Newton(func,step,Variables,nM,CauchyPointUpdated,NewtonPoint) ;  % remember also to change below in any other calls to Cauchy2Newton!

    [gammaminCN,rCN,exitflag,output]=fminbnd(FunCN,0,1,options) ;
    [stop,Outs]=outfunFminbnd();

    xVector=Outs.x ;  fVector=Outs.f ;
    I=~isnan(xVector)  ; xVector=xVector(I) ; fVector=fVector(I) ;

    if rCN < rmin

        NoReduction=false;   
   
        [r2Test,rForce,rWork,DuC2N,DvC2N,DhC2N,DlC2N]=Cauchy2Newton(func,gammaminCN,Variables,nM,CauchyPointUpdated,NewtonPoint) ;
        
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

        
        normCN=norm([du;dv;dh;dl]) ;

    end


    if CtrlVar.InfoLevelNonLinIt >= 1000
        
        
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
        yline(r0,"--k","$r_0$",interpreter="latex",LabelHorizontalAlignment="center")
        yline(rminNewton,"-.k","$\min r_N$",interpreter="latex",LabelHorizontalAlignment="right")
        yline(r2MD,"-.k","$\min r_C$",interpreter="latex",LabelHorizontalAlignment="left")
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
        funcCauchyD=@(gamma) func(gamma,dJdu,dJdv,dJdl) ;
    elseif Variables=="-uvhl-"
        funcCauchyD=@(gamma) func(gamma,dJdu,dJdv,dJdh,dJdl) ;
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

    gVector=linspace(0,2*gammaCauchyD,nPoints)' ;
    gVector(2)=gammaCauchyD/1000;

    QVector=gVector*0+nan;
    rVector=gVector*0+nan;
    for I=1:numel(gVector)

        QVector(I)=Q(gVector(I));
        rVector(I)=funcCauchyD(gVector(I)) ;
    end

    fig=FindOrCreateFigure("Q and r steepest Descent") ; clf(fig) ;

    hold on
    plot(gVector,rVector,"o-b")
    plot(gVector,QVector,"*-r")
    legend("r^2","Quad approximation in steepest Descent direction")
    plot(gammaminCauchyD,rminCauchyD,'o',MarkerFaceColor="b",MarkerSize=10)

  %% Quad approximation in the direction of M-modified steepest Descent
    % setting direction the be the that of M-modified steepest Descent
    Q=@(g) (R-g*(H*sM))'*(R-g*(H*sM))/Normalisation ;

    

    gVector=linspace(0,2*gammaminCauchyM,nPoints)' ;
    gVector(2)=gammaminCauchyM/1000;

    QVector=gVector*0+nan;
    rVector=gVector*0+nan;
    for I=1:numel(gVector)

        QVector(I)=Q(gVector(I));
        rVector(I)=funcCauchyM(gVector(I)) ;

    end

    rQM=Q(1)/funcCauchyM(1) ;

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
        r2(I)=Cauchy2Newton(func,gVector(I),"-uvhl-",nM,CauchyPointUpdated,NewtonPoint) ;
    end


    options = optimset('Display','iter','TolX',0.02,'OutputFcn',@outfunFminbnd);


    FunCN=@(step) Cauchy2Newton(func,step,"-uvhl-",nM,CauchyMPoint,NewtonPoint) ;
    [gammaminCN,rCN,exitflag,output]=fminbnd(FunCN,0,1,options) ;
    [stop,Outs]=outfunFminbnd();
    xVector=Outs.x ;  fVector=Outs.f ;
    I=~isnan(xVector)  ; xVector=xVector(I) ; fVector=fVector(I) ;

    if rCN < rmin

        [r2Test,rForce,rWork,DuC2N,DvC2N,DhC2N,DlC2N]=Cauchy2Newton(func,gammaminCN,"-uvhl-",nM,CauchyPointUpdated,NewtonPoint) ;
        
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
    r2MD=funcCauchyM(gammaminCauchyM) ;
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

rRatioMin=0.99999 ;
if NoReduction || rmin/r0 >  rRatioMin
    BackTrackInfo.Converged = false ;
    CtrlVar.InfoLevelNonLinIt = 10 ;
else
    BackTrackInfo.Converged = true ;
end

%% Summary

if CtrlVar.InfoLevelNonLinIt >= 10
    fprintf(" [---------- rLineminUa: \n")
    fprintf("\t r0=%-13.7g \t r1/r0=%-13.7g \t rNewton/r0=%-13.7g \t rminCauchyM/r0=%-13.7g \t rDescent/r0=%-13.7g \t rCN/r0=%-13.7g \n",r0,r1/r0,rminNewton/r0,rminCauchyM/r0,rminCauchyD/r0,rCN/r0)
    fprintf("\t g0=%-13.8g \t    g1=%-13.7g \t    gNewton=%-13.7g \t             gM=%-13.7g \t    gDescent=%-13.7g \t    gCM=%-13.7g \n",0,1,gammaminNewton,gammaminCauchyM/gammaCauchyM,gammaminCauchyD,gammaminCN)
    fprintf("\t      \t      \t     \t \t               \t normNewton=%-13.7g \t    normCauchyM=%-13.7g \t normCauchyD=%-13.7g \t normCN=%-13.7g \n ",...
        normNewton/normNewtonStep,normCauchyM/normNewtonStep,normCauchyD/normNewtonStep,normCN/normNewtonStep)
    if BestMethod=="Cauchy2Newton"
        fprintf("\t \t \t \t minC2N/rCauchy=%f \t minC2N/minN=%f   \n",rCN/rminCauchyM,rCN/rminNewton)
    end
    fprintf("\t =================>  Best method is %s  with rmin/r0=%g    <=====================\n",BestMethod,rmin/r0)
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


















