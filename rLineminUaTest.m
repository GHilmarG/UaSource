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
gammaminNewton=nan ; gammaminDescent=nan ; gammaminCauchy=nan ;  gammaminCauchyM=nan; gammaminCauchyM=nan ;

%% these will be the returns if no min is found
gammamin=0 ; rmin=r0 ;   du=du0*0 ; dv=dv0*0 ; dh=dh0*0 ; dl=dl0*0 ;
rForce=r0 ; rWork=nan ; D2=nan ;

%%


BestMethod="" ;

% CtrlVar.rLineMinUa="-Newton-Steepest Descent-Steepest Descent Mass-Cauchy-" ;
% CtrlVar.rLineMinUa="-Newton-Steepest Descent Mass-" ;

CtrlVar.rLineMinUa="-Auto-" ;  % First do the Newton step, ie H n = -g and evaluated r^2(x+n)/r(x)
                               % If reduction not sufficient, do Newton-direction backtracking, ie min r^2(x+ gamma n)
                               % If best backtracking step too small, ie gammaminNewton < gammaMinNewtonAccepted, do the M-Cauchy step, is solve M c = -g and evaluate r^2(x+c)/r(x) 
                               % If full Cauchy step not OK, then do Cauchy backtrackiing 
                               % If full Cauchy step OK, do minimisation along c to n direction (dog-leg)


CtrlVar.rLineMinUa="-Newton-Steepest Descent-Steepest Descent Mass-Plot Quad Approximations-" ;

CtrlVar.rLineMinUa="-Newton Step-Cauchy M-step-Steepest Descent-Plot Quad Approximations-" ;


CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ;

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
    NewtonSlope0=-2*r0 ;

    CtrlVar.uvMinimisationQuantity="Force Residuals" ;
    CtrlVar.BacktrackingGammaMin=0.001 ; 
    CtrlVar.BacktracFigName="Line Search in Newton Direction" ;
    [gammaminNewton,rminNewton,BackTrackInfo]=BackTracking(NewtonSlope0,1,r0,r1,rNewtonFunc,CtrlVar);

    if BackTrackInfo.Converged

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

        if  contains(CtrlVar.rLineMinUa,"-Auto-") &&  gammaminNewton < 0.02

            fprintf(' rLineminUa: Newton step smaller that 0.02  \n ') ;
            fprintf(' rLineminUa: Will now try using Cauchy step Steepest Descent were the mass matrix replaces the Hessian \n ') ;
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

if contains(CtrlVar.rLineMinUa,"-Cauchy M-step-") ||  contains(CtrlVar.rLineMinUa,"-Steepest Descent-") ||  contains(CtrlVar.rLineMinUa,"-Cauchy-")

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

%% M-Steepest Descent, effectivily as if using the Mass matrix as the Hessian
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

    CauchyMSlope0=(-2*R'*H*s)/Normalisation;

    gammaminCauchyM=(s'*H'*R+R'*H*s)/(2*(H*s)'*(H*s)) ;
    CauchyMPoint=gammaminCauchyM*s ;
     rD0=r0 ; 


    if CauchyMSlope0 > 0
        CauchyMSlope0=-CauchyMSlope0;
        gammaminCauchyM = -0.1 *rD0/CauchyMSlope0 ;  % initial step size
    end

    rCauchyM=rMassFunc(gammaminCauchyM);
    CtrlVar.uvMinimisationQuantity="Force Residuals" ;  CtrlVar.BacktracFigName="Line Search in Cauchy M direction" ;
    
    CtrlVar.LineSearchAllowedToUseExtrapolation=true;
    CtrlVar.NewtonAcceptRatio=0.999; % use the Armijo's condition based on slope0

    [gammaminCauchyM,rminCauchyM,BackTrackInfoCauchyM]=BackTracking(CauchyMSlope0,gammaminCauchyM,r0,rCauchyM,rMassFunc,CtrlVar);


    if BackTrackInfoCauchyM.Converged
        if  contains(CtrlVar.rLineMinUa,"-Auto-")
            if gammaminCauchyM > gammaminCauchyM / 2   % if minimum was found close to the Cauchy point, try Cauch M to Newton direction as well
                CtrlVar.rLineMinUa="-Auto-Cauchy M to Newton-" ;
            end
        end

        if rminCauchyM < rmin     % if this rmin is the smallest so far, updated solution

            du=gammaminCauchyM*Du0;
            dv=gammaminCauchyM*Dv0;
            dh=gammaminCauchyM*Dh0;
            dl=gammaminCauchyM*Dl0;
            gammamin=gammaminCauchyM;
            rmin=rminCauchyM;
            BackTrackInfo=BackTrackInfoCauchyM;
            BackTrackInfo.Direction="MD" ;
            BestMethod="Mass Steepest Descent" ;
            [rTest,~,~,rForce,rWork,D2]=rMassFunc(gammamin);
            gammamin=gammamin/gammaminCauchyM ;  % on return, normalize this with the min of the quad model

        end
    else

        fprintf(' rLineminUa: backtracking step in M-modified steepest-descent direction did not converge \n ') ;
        if  contains(CtrlVar.rLineMinUa,"-Auto-")

            fprintf(' rLineminUa: Will now try using steepest descent  \n ') ;
            CtrlVar.rLineMinUa="-Auto-Steepest Descent-";
            % If Newton step is not accepted, try Steepest Descent.

        end
    end
end



%% Steepest Descent, effectivily as if using the Unit matrix as the Hessian
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

    gammaminCauchy=0.5*R'*(H'+H)*R/(R'*(H'*H)*R) ;
    b = -0.1 *r0/slope0Descent ;  % initial step size
    b=gammaminCauchy ;

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
            BackTrackInfo.Direction="SD" ;
            BestMethod="Steepest Descent" ;
            [rTest,~,~,rForce,rWork,D2]=rDescentFunc(gammamin);
            gammamin=gammamin/gammaminCauchy ; % on return, normalize this with the min of the quad model
        end
    else

        fprintf(' rLineminUa: backtracking step in steepest-descent direction did not converge \n ') ;

    end
end

%% Dogleg, going from the M-Cauchy point towards the Newton point
if contains(CtrlVar.rLineMinUa,"-Cauchy M to Newton-")


    options = optimset('Display','iter','TolX',0.02,'OutputFcn',@outfunFminbnd);
    FunCN=@(step) Cauchy2Newton(func,step,"-uvhl-",nM,CauchyMPoint,NewtonPoint) ;
    [gammaminCN,rCN,exitflag,output]=fminbnd(FunCN,0,1,options) ;
    [stop,Outs]=outfunFminbnd();

    xVector=Outs.x ;  fVector=Outs.f ;
    I=~isnan(xVector)  ; xVector=xVector(I) ; fVector=fVector(I) ;

    if rCN < rmin

        [r2Test,rForce,rWork,DuC2N,DvC2N,DhC2N,DlC2N]=Cauchy2Newton(func,gammaminCN,"-uvhl-",nM,CauchyMPoint,NewtonPoint) ;
        
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


    if CtrlVar.InfoLevelBackTrack >= 1000
        r2Newton=rNewtonFunc(1) ;
        r2M=rMassFunc(gammaminCauchyM) ;
        figCN=FindOrCreateFigure(" Cauchy to Newton ") ; clf(figCN) ;

        [~,I]=sort(xVector)  ; xVector=xVector(I) ; fVector=fVector(I) ;

        plot(xVector,fVector,'o-r') ;
        hold on ;
        plot(0,r2M,"*m")
        plot(1,r2Newton,"*b")
        text(0,r2M,"   C(M)",HorizontalAlignment="left")
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
        dx=[du0;dv0;dl0] ;   % Newton direction
        s=[Du0;Dv0;Dl0];     % M direction
        R=[dJdu;dJdv;dJdl];  % Steepest descent
    elseif Variables=="-uvhl-"
        dx=[du0;dv0;dh0;dl0] ;    % Newton direction
        s=[Du0;Dv0;Dh0;Dl0];      % M direction
        R=[dJdu;dJdv;dJdh;dJdl];  % Steepest descent

    end

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

    %% Quad approximation in the direction of steepest descent
    dx=R ;  % setting direction the be the that of steepest descent
    Q=@(g) (R-g*(H*dx))'*(R-g*(H*dx))/Normalisation ;

    gammaminCauchy=0.5*R'*(H+H')*R/(R'*(H'*H)*R) ;

    if gammaminCauchy<0
        gammaminCauchy=-gammaminCauchy;
    end
    gVector=linspace(0,2*gammaminCauchy,nPoints)' ;
    gVector(2)=gammaminCauchy/1000;

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
    plot(gammaminDescent,rminDescent,'o',MarkerFaceColor="b",MarkerSize=10)

  %% Quad approximation in the direction of M-modified steepest descent
    dx=s ;                % setting direction the be the that of M-modified steepest descent
    Q=@(g) (R-g*(H*dx))'*(R-g*(H*dx))/Normalisation ;

    % CauchyMSlope0=(-2*R'*H*s)/Normalisation;
    % gammaminCauchyM=(s'*H'*R+R'*H*s)/(2*(H*s)'*(H*s)) ;
    % if gammaminCauchyM<0
    %    gammaminCauchyM=-gammaminCauchyM;
    % end

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
    fig=FindOrCreateFigure("Q and r M-modified steepest descent") ; clf(fig) ;

    hold on
    plot(gVector,rVector,"o-b")
    plot(gVector,QVector,"*-r")
    
    plot(gammaminCauchyM,rminCauchyM,'o',MarkerFaceColor="b",MarkerSize=10)
    plot(gammaminCauchyM,rCauchyM,'o',MarkerFaceColor="r",MarkerSize=5)
    legend("r^2","Quad approximation in M-modified steepest descent direction","Estimated min","Initial Guess (MCauchy Point)")

    %% Plot r2 along Cauchy to Newton direction


    

    gVector=linspace(0,1,nPoints)' ;
    [nM,mM]=size(M);
    r2=nan(nPoints,1) ;

    for I=1:nPoints
        r2(I)=Cauchy2Newton(func,gVector(I),"-uvhl-",nM,CauchyMPoint,NewtonPoint) ;
    end

    options = optimset('Display','iter','TolX',0.02,'OutputFcn',@outfunFminbnd);
    FunCN=@(step) Cauchy2Newton(func,step,"-uvhl-",nM,CauchyMPoint,NewtonPoint) ;
    [gammaminCN,rCN,exitflag,output]=fminbnd(FunCN,0,1,options) ;
    [stop,Outs]=outfunFminbnd();
    xVector=Outs.x ;  fVector=Outs.f ;
    I=~isnan(xVector)  ; xVector=xVector(I) ; fVector=fVector(I) ;

    if rCN < rmin

        [r2Test,rForce,rWork,DuC2N,DvC2N,DhC2N,DlC2N]=Cauchy2Newton(func,gammaminCN,"-uvhl-",nM,CauchyMPoint,NewtonPoint) ;
        
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
    r2M=rMassFunc(gammaminCauchyM) ;
    figCN=FindOrCreateFigure(" Cauchy to Newton ") ; clf(figCN) ;
    
    
    xVector=[xVector(:);gVector(:)] ; fVector=[fVector(:);r2(:)] ; 
    [~,I]=sort(xVector)  ; xVector=xVector(I) ; fVector=fVector(I) ; 

    plot(xVector,fVector,'o-r') ; 
    hold on ; 
    plot(0,r2M,"*m")
    plot(1,r2Newton,"*b")
    text(0,r2M,"   C(M)",HorizontalAlignment="left")
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

%% Summary
fprintf(" [---------- rLineminUa: \n")
fprintf("\t r0=%-13.7g \t r1/r0=%-13.7g \t rNewton/r0=%-13.7g \t rM/r0=%-13.7g \t rDescent/r0=%-13.7g \t rCauchy/r0=%-13.7g \t rCM/r0=%-13.7g \n",r0,r1/r0,rminNewton/r0,rminCauchyM/r0,rminDescent/r0,rminCauchy/r0,rCN/r0)
fprintf("\t g0=%-13.8g \t    g1=%-13.7g \t    gNewton=%-13.7g \t    gM=%-13.7g \t    gDescent=%-13.7g \t    gCauchy=%-13.7g \n",0,1,gammaminNewton,gammaminCauchyM/gammaminCauchyM,gammaminDescent,gammaminCauchy)
fprintf("\t =================>  Best method is %s     <=====================\n",BestMethod)
fprintf(" -------------------------] \n")

%%


return

end
 

 % 
 % [---------- rLineminUa: 
	%  r0=3.041337e-06  	 r1/r0=1358.289      	 rNewton/r0=0.9940949     	 rM/r0=0.7306787     	 rDescent/r0=0.8374294     	 rCauchy/r0=NaN           	 rCM/r0=0.7052111     
	%  g0=0             	    g1=1             	    gNewton=0.02910391    	    gM=1             	    gDescent=1.90946e-11   	    gCauchy=1.90946e-11   
	%  =================>  Best method is Cauchy2Newton     <=====================
 % 


















