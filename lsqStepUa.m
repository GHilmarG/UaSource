



function [R2min,dx,dlambda,gammamin,Slope0,BackTrackInfo,gammaEst,exitflag]=lsqStepUa(CtrlVar,fun,x0,lambda0,K0,R0,L,c,R20)


%
%   [H0 L' ]  [dx]  = - [g]
%   [L  0  ]  [dl]      [h]
%
%
%


nargoutchk(8,8)
narginchk(9,9)

exitflag=0 ;

isLSQ=CtrlVar.lsqUa.isLSQ ;
CostMeasure=CtrlVar.lsqUa.CostMeasure;
Step=CtrlVar.lsqUa.Step;
nx=numel(x) ;

if ~isempty(L)
    LTlambda=L'*lambda0 ;
    h0 =- (L*x0-c);
else
    LTlambda=0;
    h0=[];
end


if isLSQ
    KK0=K0'*K0;
    
    if CtrlVar.lsqUa.Step=="-Newton-"
        H0=2*KK0;
    elseif CtrlVar.lsqUa.Step=="-Cauchy-"
        H0=speye(nx) ;
    else
        error("what step?")
    end

    g0 =- (2*K0'*R0 + LTlambda) ;
else
    
    
    if CtrlVar.lsqUa.Step=="-Newton-"
        H0=K0;
    elseif CtrlVar.lsqUa.Step=="-Cauchy-"
        H0=speye(nx) ;
    else
        error("what step?")
    end

    KK0=K0'*K0;
    g0 =- (R0 + LTlambda) ;

end





CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical=issymmetric(K0) ;

if CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical
    [dx,dlambda]=solveKApeSymmetric(H0,L,g0,h0,x0,lambda0,CtrlVar);
else
    [dx,dlambda]=solveKApe(H0,L,g0,h0,x0,lambda0,CtrlVar);
end


if CostMeasure=="R2"

    J0=full(R0'*R0) ;
    Slope0=2*R0'*(K0*dx) ;
    gammaEst=-(R0'*K0*dx)/(dx'*(KK0)*dx) ;

elseif CostMeasure=="r2"
    
    
    d=[g0;h0] ;
    J0=full(d'*d);
    Hd=[H0*dx+L'*dlambda,L*dx];   % this should be equal to -d
    Slope0=2*d'*Hd;               % this should be equal to -2*d'*d 
    gammaEst=d'*Hd/(Hd'*Hd) ;     % I have the minus in the solve

end



if Slope0 > 0

    %fprintf("lsqStepUa: Exiting because slope at origin in line search positive (Slope=%g) \n",Slope0)
    R2min=nan;
    gammamin=nan ; BackTrackInfo=[];
    gammaEst=nan ;
    exitflag=1 ;
    return
end




CtrlVar.BacktrackingGammaMin=gammaEst*CtrlVar.BacktrackStepRatio ;

funcBackTrack=@(gamma) Jlsqfunc(CtrlVar,gamma,dx,dlambda,fun,L,c,x0,lambda0) ;


R2=nan;

CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ;  CtrlVar.doplots=1 ;

%CtrlVar.NewtonAcceptRatio=0.001;


[gammamin,R2min,BackTrackInfo]=BackTracking(Slope0,gammaEst,R20,R2,funcBackTrack,CtrlVar);

% dx=gammamin*dx  ;
% dlambda=gammamin*dlambda ;
% x=x0+dx ;
% lambda=lambda0+dlambda ;
% R2=R2min ;


end