



function [Jmin,dx,dlambda,gammamin,Slope0,BackTrackInfo,gammaEst,exitflag]=lsqStepUa(CtrlVar,fun,x0,lambda0,K0,R0,L,c)


%
%   [H0 L' ]  [dx]  = - [g]
%   [L  0  ]  [dl]      [h]
%
%
%


nargoutchk(8,8)
narginchk(8,8)

exitflag=0 ;

isLSQ=CtrlVar.lsqUa.isLSQ ;
CostMeasure=CtrlVar.lsqUa.CostMeasure;
Step=CtrlVar.lsqUa.Step;
nx=numel(x0) ;

if ~isempty(L)
    LTlambda=L'*lambda0 ;
    h0 =- (L*x0-c);
else
    LTlambda=0;
    h0=[];
end


if isLSQ
    if CtrlVar.lsqUa.Step=="-Newton-"
        H0=2*(K0'*K0);
    elseif CtrlVar.lsqUa.Step=="-Cauchy-"
        H0=speye(nx) ;
    else
        error("what step?")
    end
    g0 =- (2*K0'*R0 + LTlambda) ;
else
    if Step=="-Newton-"
        H0=K0;
    elseif Step=="-Cauchy-"
        H0=speye(nx) ;
    else
        error("what step?")
    end
    g0 =- (R0 + LTlambda) ;
end

CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical=issymmetric(K0) ;

if CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical
    [dx,dlambda]=solveKApeSymmetric(H0,L,g0,h0,x0,lambda0,CtrlVar);
else
    [dx,dlambda]=solveKApe(H0,L,g0,h0,x0,lambda0,CtrlVar);
end

if CostMeasure=="R2"

    % Here gammaEst should always be equal to 1 for Newton only.
    % Also for Cauchy in an unconstrained case 

    J0=full(R0'*R0) ;
    K0dx=K0*dx ;
    Slope0=full(2*R0'*K0dx) ;
    gammaEst=-full((R0'*K0dx)/(K0dx'*K0dx));

elseif CostMeasure=="r2"
    
    %
    %  Here gammaEst should always be equal to 1, both for Newton and Cauchy
    %

    d=[g0;h0] ;
    J0=full(d'*d);
    if ~isempty(L)
        Hd=[H0*dx+L'*dlambda; L*dx];   % this should be equal to -d
    else
        Hd=H0*dx ;
    end
    Slope0=-full(2*d'*Hd);              % this should be equal to -2*d'*d
    gammaEst=full(d'*Hd/(Hd'*Hd)) ;     % I have the minus in the solve

end

s=-2*J0; 
% fprintf("Slope0=%g \t -2J0=%g \t Slope0/(-2J0)=%g \n ",Slope0,s,Slope0/s)



if Slope0 > 0

    %fprintf("lsqStepUa: Exiting because slope at origin in line search positive (Slope=%g) \n",Slope0)
    Jmin=nan;
    gammamin=nan ; BackTrackInfo=[];
    gammaEst=nan ;
    exitflag=1 ;
    return
end


CtrlVar.BacktrackingGammaMin=gammaEst*CtrlVar.BacktrackStepRatio ;

funcBackTrack=@(gamma) Jlsqfunc(CtrlVar,gamma,dx,dlambda,fun,L,c,x0,lambda0,K0) ;

J=nan;

% CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ;  CtrlVar.doplots=1 ;

%CtrlVar.NewtonAcceptRatio=0.001;

[gammamin,Jmin,BackTrackInfo]=BackTracking(Slope0,gammaEst,J0,J,funcBackTrack,CtrlVar);


Jmin=full(Jmin) ;
gammamin=full(gammamin) ;
Slope0=full(Slope0); 
gammaEst=full(gammaEst);



end