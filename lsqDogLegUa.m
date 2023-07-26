function [x,lambda,R2,g2,residual,g,h,output] = lsqDogLegUa(CtrlVar,fun,x,lambda,L,c)




if isempty(CtrlVar) || ~isstruct(CtrlVar)

    ItMax=5;
    gTol=1e-20;
    dR2Tol=1e-3;
    dxTol=1e-20;

    isLSQ=true;
    Normalize=false ;
    SaveIterate=true;
    LevenbergMarquardt="auto" ; % "fixed"
    LMlambda=1 ;
    ScaleProblem=false;
    LMlambdaUpdateMethod=1;


else

    ItMax=CtrlVar.lsqUa.ItMax ;
    gTol=CtrlVar.lsqUa.gTol ;
    dR2Tol=CtrlVar.lsqUa.dR2Tol ;
    dxTol=CtrlVar.lsqUa.dxTol ;

    isLSQ=CtrlVar.lsqUa.isLSQ ;
    Normalize=CtrlVar.lsqUa.Normalize;
    LevenbergMarquardt=CtrlVar.lsqUa.LevenbergMarquardt;
    LMlambda=CtrlVar.lsqUa.LMlambda0 ;
    ScaleProblem=CtrlVar.lsqUa.ScaleProblem;

    LMlambdaUpdateMethod=CtrlVar.lsqUa.LMlambdaUpdateMethod ;

    SaveIterate=CtrlVar.lsqUa.SaveIterate;

end

R2Array=nan(ItMax+1,1) ;
g2Array=nan(ItMax+1,1) ;
dR2=[inf ; inf ] ; % stores the changes in R2=R'*R  over last two iterations

nx=numel(x);

%% If contraints provided, make iterate feasable

if ~isempty(L) && ~isempty(c)   % if the user has not provided an initial estimate of lambda, but specifies constraints, set lambda=0

    BCres=norm(L*x-c);
    if BCres>1e-6   % make feasable
        x=L\c ;
    end
    if isempty(lambda)  || isnan(lambda)
        lambda=c*0;
    end
end

if ~isempty(L)
    LTlambda=L'*lambda ;
else
    LTlambda=0;
end
%%

%% Normalisation
if Normalize
    x0=x*0 ;
    R=fun(x0) ;

    g =- (R + LTlambda) ;
    Normalisation=full(g'*g);
else
    Normalisation=1;

end

% Evaluate cost function, don't solve system
[R,K]=fun(x) ;

if isLSQ
    g =- (2*K'*R + LTlambda) ;
else
    g =- (R + LTlambda) ;
end

if ~isempty(L)
    h =- (L*x-c);
else
    h=[];
end

R2=full(R'*R);
g2 = full(g'*g)/Normalisation;


g2Array(1)=g2;
R2Array(1)=R2;
if SaveIterate
    xVector=nan(numel(x),100);
    xVector(:,1)=x(:) ;
else
    xVector=[];
end
iteration=0 ;

fprintf("\n\t Start lsqUa: \t  g=%g \t         r=%g \n \n",g2,R2)


%%


while iteration <= ItMax

    iteration=iteration+1 ;

    K0=K ; R0=R; x0=x ; lambda0=lambda ; h0=h ; g0=g;
    R20=R2;  g20=g2 ; 
    R2=nan; g2=nan ;

    if isLSQ
        KK0=K0'*K0;
        H0=2*KK0;

    else
        H0=K0;
        KK0=K0'*K0;
    end

    [dxN,dlambdaN]=solveKApe(H0,L,g0,h0,x0,lambda0,CtrlVar);


    funcBackTrack=@(gamma) R2func(gamma,dxN,dlambdaN,fun,x0,lambda0) ;


    Slope0=2*R0'*K0*dxN ;
    gammaEst=-R0'*K0*dxN/(dxN'*(KK0)*dxN) ;

    if Slope0 > 0
        fprintf("lsqUa: Exiting iteration because slope at origin in line search positive (Slope=%g) \n",Slope0)
        R2=R20 ; g2=g20 ; 
        break

    end



    %CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ; CtrlVar.NewtonAcceptRatio=0.001; CtrlVar.doplots=1 ;
    [gammaminNewton,R2minNewton,BackTrackInfo]=BackTracking(Slope0,gammaEst,R20,R2,funcBackTrack,CtrlVar);

    R2=R2minNewton ;

    dx=gammaminNewton*dxN  ;
    dlambda=gammaminNewton*dlambdaN ;
    x=x+dx ;
    lambda=lambda+dlambda ;

    if isLSQ
        Q=2*R0'*K0*dx+dx'*KK0*dx ;  % Quad approximation, based on unperturbed H
    else
        Q=R0'*dx+dx'*H0*dx/2 ;
    end


    [R,K]=fun(x) ;
    R2=full(R'*R);

    if ~isempty(L)
        LTlambda=L'*lambda ;
        h =- (L*x-c);
    else
        LTlambda=0;
        h=[];
    end

    if isLSQ
        g =- (2*K'*R + LTlambda) ;
    else
        g =- (R + LTlambda) ;
    end

    g2=full(g'*g)/Normalisation ;
    
    
   

    rho=(R2-R20)/Q;      % Actual reduction / Modelled Reduction
    g2Ratio=g2/g20 ;
    R2Ratio=R2/R20 ;
    dR2=[abs(R2-R20); dR2(1)] ;


    dxNorm=norm(dx);
    dlambdaNorm=norm(dlambda);
    BCsNorm=norm(h) ;


    g2Array(iteration+1)=g2;
    R2Array(iteration+1)=R2;
    if SaveIterate
        xVector(:,iteration+1)=x(:) ;
    end


    fprintf("lsqUa: \t it=%2i  \t     |R|^2=%-13g \t     R1/R0=%-13g \t gamma=%-13g \t |g|^2=%-13g \t |dx|=%-13g \t |dl|=%-13g \t |BCs|=%-13g \t dr/Q=%-5f \t slope0 =%g \n",iteration,R2,R2Ratio,gammaminNewton,g2,dxNorm,dlambdaNorm,BCsNorm,rho,Slope0)


    if g2 < gTol
        fprintf("lsqUa: Exiting iteration because |g|^2 within set tolerance of %g \n",gTol)
        break
    end

    if dxNorm < dxTol
        fprintf("lsqUa: Exiting iteration because change in x within the set tolerance of %g \n",dxTol)
        break
    end


    maxdR2=max(dR2);
    if maxdR2 < dR2Tol
        fprintf("lsqUa: Exiting iteration because max change in |R|^2 over last two iterations (%g) less than the set tolerance of %g \n",maxdR2,dR2Tol)
        break
    end

    if iteration >= ItMax
        fprintf("lsqUa: Exiting iteration because number of iterations has reached the set maximum of %i \n",ItMax)
        break

    end


end

fprintf("\n\t Exit lsqUa: \t  g=%g \t         r=%g \n \n",g2,R2)


residual=R ;
output.g2Array=g2Array;
output.R2Array=R2Array;
output.xVector=xVector;
output.nIt=iteration;

end

function  R2=R2func(gamma,dx,dl,fun,x0,l0)

x=x0+gamma*dx;
l=l0+gamma*dl ;

R=fun(x) ;
R2=R'*R;



end




