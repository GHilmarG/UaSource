



function [Jmin,dx,dlambda,gammamin,Slope0,BackTrackInfo,gammaEst,exitflag]=lsqStepUa(CtrlVar,fun,x0,lambda0,K0,R0,L,c)

%%
%
% Newton:
%
%   [K0 L' ]  [dx]  = - [R0 - L' lambda0]
%   [L  0  ]  [dl]      [L x0 - c        ]
%
%
%
% Newton lsq
%
%   [2 K0'K0   L' ]  [dx]  = - [2 K0' R0 + L' lambda0 ]
%   [    L     0  ]  [dl]      [      L x0 - c        ]
%
%
%
% Cauchy:
%
%   [    I     L' ]  [dx]  = - [2 K0' R0 + L' lambda0 ]
%   [    L     0  ]  [dl]      [      L x0 - c        ]
%
% search direction:  d=[dx ; dlambda]
%
% Minimizing R^2 where I know
%
%   H = \nabla R
%
% in some direction d
%
% lsq:
%
%  Q=2*R'*K*dx+dx'*KK*dx ; 
%  Q(gamma)=2 gamma (R' K)' dx + gamma^2 dx' KK dx 
%  dQ/dgamma= 2 (R' K)' dx + 2 gamma dx' KK dx
%
%  Slope0= 2 R' K dx
%  gammaMinQ = - R' K dx /( (dx' KK dx )
%
%
% ~lsq
%    Q=R'*dx+dx'*H*dx/2 ;
%
%%


nargoutchk(8,8)
narginchk(8,8)

exitflag=0 ;

isLSQ=CtrlVar.lsqUa.isLSQ ;
CostMeasure=CtrlVar.lsqUa.CostMeasure;

nx=numel(x0) ;

if ~isempty(L)
    LTlambda=L'*lambda0 ;
    h0 =- (L*x0-c);
else
    LTlambda=0;
    h0=[];
end


if isLSQ
    H0=2*(K0'*K0);
    g0 =- (2*K0'*R0 + LTlambda) ;
else

    H0=K0;
    g0 =- (R0 + LTlambda) ;
end


if CtrlVar.lsqUa.Step=="-Newton-"

    CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical=issymmetric(H0) ; % this should not really be needed as this will always be a symmetrical pos def matrix

    if CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical
        [dx,dlambda]=solveKApeSymmetric(H0,L,g0,h0,x0,lambda0,CtrlVar);
    else
        [dx,dlambda]=solveKApe(H0,L,g0,h0,x0,lambda0,CtrlVar);
    end

elseif CtrlVar.lsqUa.Step=="-Cauchy-"

    % This is a simple system. I should be able to figure out the solution
    % without calling the solver
    g0= -(2*K0'*R0 + LTlambda) ;
    HCauchy=speye(nx) ;
    CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical=issymmetric(HCauchy) ;

    if CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical
        [dx,dlambda]=solveKApeSymmetric(HCauchy,L,g0,h0,x0,lambda0,CtrlVar);
    else
        [dx,dlambda]=solveKApe(HCauchy,L,g0,h0,x0,lambda0,CtrlVar);
    end
    % Cauchy not working, presumably a sign error in the direction

else
    error("case not found")
end

if CostMeasure=="R2"

    % Here gammaEst should always be equal to 1 for Newton only.
    % Also for Cauchy in an unconstrained case

    J0=full(R0'*R0) ;
    K0dx=K0*dx ;
    Slope0=full(2*R0'*K0dx) ;
    gammaEst=-full((R0'*K0dx)/(K0dx'*K0dx));

elseif CostMeasure=="r2"

    r=[g0 ; h0 ] ;
    J0=full(r'*r) ;


    if ~isempty(L)
        Hd=[H0*dx+L'*dlambda; L*dx];   % this should be equal to -d
    else
        Hd=H0*dx ;
    end

    Slope0=-full(2*r'*Hd);              % this should be equal to -2*r'*r
    gammaEst=full(r'*Hd/(Hd'*Hd)) ;     % I have the minus in the solve

    % Cauchy: Q=g x + 0.5 x H x
    % minimize in the direction -g with respect to alpha
    %
    % J=-alpha g' g + alpha^2 0.5 g H g
    %
    % dJ/dalpha= -g' g + alpha g H g - > alpha = g'g/(gHg)
    H=[H0 L' ; L sparse(size(L,1),size(L,1))] ;
    gammaEstCauchy=r'*r/(r'*H*r) ;

    if CtrlVar.lsqUa.Step=="-Cauchy-"
        gammaEst=gammaEstCauchy;
    end

    if gammaEst <0
        gammaEst=-0.9/Slope0;
    end
end



% s=-2*J0;
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

funcBackTrack=@(gamma) Jlsqfunc(CtrlVar,gamma,dx,dlambda,fun,L,c,x0,lambda0) ;

J=nan;
J0=funcBackTrack(0) ; % is this same as J0 above?

CtrlVar.InfoLevelBackTrack=10000;  CtrlVar.InfoLevelNonLinIt=10 ;  CtrlVar.doplots=1 ;

%CtrlVar.NewtonAcceptRatio=0.001;

[gammamin,Jmin,BackTrackInfo]=BackTracking(Slope0,gammaEst,J0,J,funcBackTrack,CtrlVar);


Jmin=full(Jmin) ;
gammamin=full(gammamin) ;
Slope0=full(Slope0);
gammaEst=full(gammaEst);



end