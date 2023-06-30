function [gammamin,rmin,du,dv,dh,dl,BackTrackInfo,rForce,rWork,D2] = rLineminUaTest(CtrlVar,UserVar,func,r0,r1,K,L,du0,dv0,dh0,dl0,dJdu,dJdv,dJdh,dJdl,Normalisation,M)

%%

%  [ K  L ]  [dx0]  = [ dJdx ]
%  [ L' 0 ]  [dl0]    [ dJdl ]
%


% func=@(gamma,Du,Dv,Dl) CalcCostFunctionNR(UserVar,RunInfo,CtrlVar,MUA,gamma,F,fext0,L,l,cuv,Du,Dv,Dl) ;
%%

% CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ;

%%

rminNewton=nan     ; rminDescent=nan     ; rminCauchy=nan     ;  rminM=nan;
gammaminNewton=nan ; gammaminDescent=nan ; gammaMinCauchy=nan ;  gammaminM=nan; gammaMinCauchyM=nan ;

%% these will be the returns if no min is found
gammamin=0 ; rmin=r0 ;   du=du0*0 ; dv=dv0*0 ; dh=dh0*0 ; dl=dl0*0 ;
rForce=r0 ; rWork=nan ; D2=nan ;

%%


% CtrlVar.rLineMinUa="-Newton-Steepest Descent-Steepest Descent Mass-Cauchy-" ;
% CtrlVar.rLineMinUa="-Newton-Steepest Descent Mass-" ;

CtrlVar.rLineMinUa="-Auto-" ;

CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ;

%%

if isempty(dh0)
    Variables="-uvl-";
else
    Variables="-uvhl-";
end


if contains(CtrlVar.rLineMinUa,"-Newton-")  || contains(CtrlVar.rLineMinUa,"-Auto-")

    % Newton

    % Newton direction is [du0;dv0;dl0] ;

    if Variables=="-uvl-"
        rNewtonFunc=@(gamma) func(gamma,du0,dv0,dl0) ;
    elseif Variables=="-uvhl-"
        rNewtonFunc=@(gamma) func(gamma,du0,dv0,dh0,dl0) ;
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
    slope0=-2*r0 ;

    CtrlVar.uvMinimisationQuantity="Force Residuals" ;
    CtrlVar.BacktrackingGammaMin=0.02 ; 
    CtrlVar.BacktracFigName="Line Search in Newton Direction" ;
    [gammaminNewton,rminNewton,BackTrackInfo]=BackTracking(slope0,1,r0,r1,rNewtonFunc,CtrlVar);

    if BackTrackInfo.Converged

        du=gammaminNewton*du0;
        dv=gammaminNewton*dv0;
        dh=gammaminNewton*dh0;
        dl=gammaminNewton*dl0;

        gammamin=gammaminNewton;
        rmin=rminNewton;
        [rTest,~,~,rForce,rWork,D2]=rNewtonFunc(gammamin);
        BackTrackInfo.Direction="Newton" ;

        if abs(rTest-rmin)>1000*eps

            error("rLineMinUa:inconsistent","r not having the value expected")

        end

        if  contains(CtrlVar.rLineMinUa,"-Auto-") &&  gammaminNewton < 0.02

            fprintf(' rLineminUa: Newton step smaller that 0.02  \n ') ;
            fprintf(' rLineminUa: Will now try using Steepest Descent were the mass matrix replaces the Hessian \n ') ;
            CtrlVar.rLineMinUa="-Auto-Steepest Descent Mass-Steepest Descent-" ;

        end

    else

        fprintf(' rLineminUa: backtracking step in Newton direction did not converge \n ') ;

        if  contains(CtrlVar.rLineMinUa,"-Auto-")

            fprintf(' rLineminUa: Will now try using Steepest Descent were the mass matrix replaces the Hessian \n ') ;
            CtrlVar.rLineMinUa="-Auto-Steepest Descent Mass-Steepest Descent-" ;
             

        end

    end



end

if contains(CtrlVar.rLineMinUa,"-Steepest Descent Mass-") ||  contains(CtrlVar.rLineMinUa,"-Steepest Descent-") ||  contains(CtrlVar.rLineMinUa,"-Cauchy-")

    if Variables=="-uvl-"
        R=[dJdu;dJdv;dJdl];
    elseif Variables=="-uvhl-"
        R=[dJdu;dJdv;dJdh;dJdl];
    else
        error(" Variables? ")
    end


    [nL,mL]=size(L);
    L0=sparse(nL,nL);
    H=[K L' ;L L0] ;

end

%% Descent/Newton, ie using the Mass matrix as Hessian



if contains(CtrlVar.rLineMinUa,"-Steepest Descent Mass-")

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
        [sol,Dl0]=solveKApeSymmetric(Mblock,L,[dJdu;dJdv],dJdl,[du0;dv0],dl0,CtrlVar);
        Du0=sol(1:nM) ; Dv0=sol(nM+1:2*nM);
        rMassFunc=@(gamma) func(gamma,Du0,Dv0,Dl0) ;
        s=[Du0;Dv0;Dl0];

    elseif Variables=="-uvhl-"


        Mblock=[M O O ; O M O ; O O M];
        [sol,Dl0]=solveKApeSymmetric(Mblock,L,[dJdu;dJdv;dJdh],dJdl,[du0;dv0;dh0],dl0,CtrlVar);
        Du0=sol(1:nM) ; Dv0=sol(nM+1:2*nM);  Dh0=sol(2*nM+1:3*nM); 
        rMassFunc=@(gamma) func(gamma,Du0,Dv0,Dh0,Dl0) ;
        s=[Du0;Dv0;Dh0;Dl0];

    end

    slopeD0=(-2*R'*H*s)/Normalisation;

    gammaMinCauchyM=(s'*H'*R+R'*H*s)/(2*(H*s)'*(H*s)) ;

     rD0=r0 ; 


    if slopeD0 > 0
        slopeD0=-slopeD0;
        gammaMinCauchyM = -0.1 *rD0/slopeD0 ;  % initial step size
    end

    %rD0=rMassFunc(0) ;

    % gamma=b ; rMb=rMassFunc(gamma);
    b=gammaMinCauchyM ; rMb=rMassFunc(b);


    CtrlVar.uvMinimisationQuantity="Force Residuals" ;  CtrlVar.BacktracFigName="Line Search in Mass Direction" ;
    CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ;
    CtrlVar.LineSearchAllowedToUseExtrapolation=true;
    CtrlVar.NewtonAcceptRatio=0.999; % use the Armijo's condition based on slope0

    [gammaminM,rminM,BackTrackInfoSteepestMass]=BackTracking(slopeD0,b,r0,rMb,rMassFunc,CtrlVar);

    if BackTrackInfoSteepestMass.Converged
        if rminM < rmin     % if this rmin is the smallest so far, updated solution
            du=gammaminM*Du0;
            dv=gammaminM*Dv0;
            dh=gammaminM*Dh0;
            dl=gammaminM*Dl0;
            gammamin=gammaminM;
            rmin=rminM;
            BackTrackInfo=BackTrackInfoSteepestMass;
            BackTrackInfo.Direction="Mass Steepest Descent" ;
            [rTest,~,~,rForce,rWork,D2]=rMassFunc(gammamin);
            gammamin=gammamin/gammaMinCauchyM ;  % on return, normalize this with the min of the quad model
        end
    else

        fprintf(' rLineminUa: backtracking step in M-modified steepest-descent direction did not converge \n ') ;

        if  contains(CtrlVar.rLineMinUa,"-Auto-")

            fprintf(' rLineminUa: Will now try using steepest descent  \n ') ;
            CtrlVar.rLineMinUa="-Auto-Steepest Descent Mass-Steepest Descent-";
            % If Newton step is not accepted, try both M-Descent and Steepest Descent.  

        end


    end

   

end

%% descent


if contains(CtrlVar.rLineMinUa,"-Steepest Descent-")

    slope0Descent=-2*R'*H*R/Normalisation ;
    % This slope estimate depends on H>0. If that is not the case
    % I simply change the sign. This will not be accurate
    if slope0Descent > 0
        slope0Descent = - slope0Descent;
    end


    if Variables=="-uvl-"
        rDescentFunc=@(gamma) func(gamma,dJdu,dJdv,dJdl) ;
    elseif Variables=="-uvhl-"
        rDescentFunc=@(gamma) func(gamma,dJdu,dJdv,dJdh,dJdl) ;
    else
        error("variables?")
    end

    gammaMinCauchy=R'*H*R/(R'*(H'*H)*R) ;
    b = -0.1 *r0/slope0Descent ;  % initial step size
    b=gammaMinCauchy ;

    gamma=b ; rb=rDescentFunc(gamma);

    % CtrlVar.InfoLevelBackTrack=1000;
    CtrlVar.BacktracFigName="Steepest Descent" ;
    CtrlVar.BacktrackingGammaMin=1e-13; CtrlVar.LineSearchAllowedToUseExtrapolation=true;
    CtrlVar.NewtonAcceptRatio=0.999; % primarirly use the Armijo's condition based on slope0

    [gammaminDescent,rminDescent,BackTrackInfoSteepest]=BackTracking(slope0Descent,b,r0,rb,rDescentFunc,CtrlVar);

    if Variables=="-uvl-"
        dNewton=[du0;dv0;dl0]; 
        dSteepest=[dJdu;dJdv;dJdl];
    elseif Variables=="-uvhl-"
        dNewton=[du0;dv0;dh0;dl0]; 
        dSteepest=[dJdu;dJdv;dJdh;dJdl];
    else
        error(" Variables? ")

    end


    p=(dNewton'*dSteepest)/(norm(dNewton)*norm(dSteepest)) ;

    fprintf("angle between Newton and Steepest descent directions = %g (deg) \n",acosd(p))

    if BackTrackInfoSteepest.Converged
        if rminDescent < rmin
            du=gammaminDescent*dJdu;
            dv=gammaminDescent*dJdv;
            dh=gammaminDescent*dJdh;
            dl=gammaminDescent*dJdl;
            gammamin=gammaminDescent;
            rmin=rminDescent;
            BackTrackInfo=BackTrackInfoSteepest;
            BackTrackInfo.Direction="Steepest Descent" ;
            [rTest,~,~,rForce,rWork,D2]=rDescentFunc(gammamin);
            gammamin=gammamin/gammaMinCauchy ; % on return, normalize this with the min of the quad model
        end
    else

        fprintf(' rLineminUa: backtracking step in steepest-descent directoin did not converge \n ') ;

    end
end




%% Cauchy
if contains(CtrlVar.rLineMinUa,"-Cauchy-")

    if Variables=="-uvl-"
        rDescentFunc=@(gamma) func(gamma,dJdu,dJdv,dJdl) ;
    elseif Variables=="-uvhl-"
        rDescentFunc=@(gamma) func(gamma,dJdu,dJdv,dJdh,dJdl) ;
    else
        error("variables?")
    end

    


    gammaMinCauchy=R'*H*R/(R'*(H'*H)*R) ;

    if gammaMinCauchy<0
        gammaMinCauchy=-gammaMinCauchy;
    end

    rminCauchy=rDescentFunc(gammaMinCauchy) ;

    if rminCauchy < rmin

        du=gammaMinCauchy*dJdu;
        dv=gammaMinCauchy*dJdv;
        dh=gammaminDescent*dJdh;
        dl=gammaMinCauchy*dJdl;
        gammamin=gammaMinCauchy;
        rmin=rminCauchy;

        BackTrackInfo.Infovector=[0 r0 ; gammamin rmin] ;
        BackTrackInfo.Converged=1;
        BackTrackInfo.iarm=0;

        BackTrackInfo.Direction="Cauchy" ;

        [rTest,~,~,rForce,rWork,D2]=rDescentFunc(gammamin);


    end


end


if contains(CtrlVar.rLineMinUa,"-Plot Quad Approximations-")

    if Variables=="-uvl-"
        rDescentFunc=@(gamma) func(gamma,dJdu,dJdv,dJdl) ;
    elseif Variables=="-uvhl-"
        rDescentFunc=@(gamma) func(gamma,dJdu,dJdv,dJdh,dJdl) ;
    else
        error("variables?")
    end

    if Variables=="-uvl-"
        % Quad approximation in Newton direction
        dx=[du0;dv0;dl0] ;   % Newton direction
        R=[dJdu;dJdv;dJdl];  % Steepest descent
    elseif Variables=="-uvhl-"
        dx=[du0;dv0;dh0;dl0] ;   % Newton direction
        R=[dJdu;dJdv;dJdh;dJdl];  % Steepest descent

    end

    Q=@(g) (R-g*(H*dx))'*(R-g*(H*dx))/Normalisation ;

    gVector=linspace(0,1.2,50)' ;
    QVector=gVector*0+nan;
    rVector=gVector*0+nan;

    for I=1:numel(gVector)

        QVector(I)=Q(gVector(I));
        rVector(I)=rNewtonFunc(gVector(I)) ;
    end

    fig=FindOrCreateFigure("Q and r (Newton)") ; clf(fig) ;

    hold on
    plot(gVector,rVector,"o-b")
    plot(gVector,QVector,"*-r")


    legend("r^2","Quad approximation in Newton direction")


    % Quad approximation in the direction of steepest descent
    dx=R ;  % setting direction the be the that of steepest descent
    Q=@(g) (R-g*(H*dx))'*(R-g*(H*dx))/Normalisation ;

    gammaMinCauchy=R'*H*R/(R'*(H'*H)*R) ;

    if gammaMinCauchy<0
        gammaMinCauchy=-gammaMinCauchy;
    end
    gVector=linspace(0,2*gammaMinCauchy,50)' ;
    gVector(2)=gammaMinCauchy/1000;

    QVector=gVector*0+nan;
    rVector=gVector*0+nan;
    for I=1:numel(gVector)

        QVector(I)=Q(gVector(I));
        rVector(I)=rDescentFunc(gVector(I)) ;
    end

    fig=FindOrCreateFigure("Q and r steepest descent") ; clf(fig) ;

    hold on
    plot(gVector,rVector,"o-b")
    plot(gVector,QVector,"*-r")
    legend("r^2","Quad approximation in steepest descent direction")

end

%% Summary
fprintf(" [---------- rLineminUa: \n")
fprintf("\t r0=%-13.7g \t r1/r0=%-13.7g \t rNewton/r0=%-13.7g \t rM/r0=%-13.7g \t rDescent/r0=%-13.7g \t rCauchy/r0=%-13.7g \n",r0,r1/r0,rminNewton/r0,rminM/r0,rminDescent/r0,rminCauchy/r0)
fprintf("\t g0=%-13.8g \t    g1=%-13.7g \t    gNewton=%-13.7g \t    gM=%-13.7g \t    gDescent=%-13.7g \t    gCauchy=%-13.7g \n",0,1,gammaminNewton,gammaminM/gammaMinCauchyM,gammaminDescent,gammaMinCauchy)
fprintf(" -------------------------] \n")

%%


return

end



















