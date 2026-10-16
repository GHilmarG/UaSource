function [rForce,rWork,D2,rBlocks]=uvhResidualMeasures(CtrlVar,MUA,frhs,grhs,cuvh,fext0,dub,dvb,dh,dl)

%% Residual measures of the uvh system
%
%   [rForce,rWork,D2,rBlocks]=uvhResidualMeasures(CtrlVar,MUA,frhs,grhs,cuvh,fext0,dub,dvb,dh,dl)
%
% (15 Oct 2026) Moved here, unchanged, from CalcCostFunctionNRuvh.m, so that these measures can also be calculated directly from an
% already assembled residual (eg at the start of a Newton-Raphson iteration, gamma=0), without a further assembly.
%
% frhs, grhs : the right-hand sides of the Newton system, frhs=-R-L'*luvh and grhs=cuvh-L*x (or frhs=-R and grhs=[] without constraints)
% dub,dvb,dh,dl : the step
% rForce : the force residual, which is the merit function of the line search
% rWork  : the work residual, rWork=D2^2, used in the convergence criteria (not as a merit function)
% D2     : the Newton decrement, D2=[frhs;grhs]'*[dub;dvb;dh;dl]
% rBlocks: the contributions of the u, v, h and l blocks to rForce

D2=[frhs;grhs]'*[dub;dvb;dh;dl];
rWork=full(D2^2);

%% Force residual
%
% Two normalisations are available, selected by
%
%    CtrlVar.uvhResidualNormalisation = "pooled"    (default, historical)
%                                     = "blockwise"
%
% "pooled" divides the whole residual by the norm of the zero-field
% reference residual fext0.  The difficulty is that the u, v and h rows of
% the uvh residual are not commensurate: the momentum rows have units of
% force while the thickness rows have units of mass, and the thickness rows
% carry factors of dt. 
%
% "blockwise" weights each block by its own scale taken from the
% corresponding block of fext0,
%
%     rForce = (|Ru|/suv)^2 + (|Rv|/suv)^2 + (|Rh|/sh)^2 + (|Rl|/sl)^2
%
% which is ||W*[frhs;grhs]||^2 for a fixed block-diagonal W.  Because W is
% fixed through a Newton solve, the identity
%
%     dr/dgamma|_0 = -2 r(0)
%
% still holds exactly, so the Newton direction remains a descent direction
% and the line search is unaffected.  Unlike max(...) over blocks, the
% weighted 2-norm is smooth, which the backtracking model requires.
%
% The optional seventh output returns the individual block residuals so that
% the caller can report which equation is the least converged.
%%

N=MUA.Nnodes;

if isfield(CtrlVar,"uvhResidualNormalisation") && CtrlVar.uvhResidualNormalisation=="blockwise"

    % block scales, fixed for the whole Newton solve because fext0 and cuvh are

    suv=norm(fext0(1:2*N));      % same normalization for both u and v components
    % Modification (8 Oct 2026): geometry-independent floor for the momentum scale, so that suv
    % does not vanish for h=0 and/or zero surface slope. The floor is the l2 norm of the nodal load
    % vector of a uniform stress density tauFloor (in units of stress, e.g. kPa), ie tauFloor*sqrt(2)*||M*1||.
    if isfield(CtrlVar,"uvhResidualNormalisationTauFloor") ; tauFloor=CtrlVar.uvhResidualNormalisationTauFloor ; else ; tauFloor=0.01 ; end
    if isfield(MUA,"M") && ~isempty(MUA.M) ; Mlump=MUA.M*ones(N,1) ; else ; Mlump=MassMatrix2D1dof(MUA)*ones(N,1) ; end
    suv=max(suv,tauFloor*sqrt(2)*norm(Mlump)) ;
    sh=norm(fext0(2*N+1:3*N));   % and a separate one for the h block 

    % a block with no reference scale is measured absolutely rather than
    % being allowed to blow up
    Small=1000*eps;
    if suv<Small ; suv=1 ; end
 
    if sh<Small ; sh=1 ; end

    ru=norm(frhs(1:N))/suv;
    rv=norm(frhs(N+1:2*N))/suv;
    rh=norm(frhs(2*N+1:3*N))/sh;

    if isempty(grhs)
        rl=0;
    else
        sl=norm(cuvh);          % constant during the Newton solve
        if sl<Small ; sl=1 ; end
        rl=norm(grhs)/sl;
    end

    rForce=full(ru^2+rv^2+rh^2+rl^2);
    rBlocks=[ru rv rh rl];
    rBlocks=rBlocks .^2; 

else

    rForce=full([frhs;grhs]'*[frhs;grhs]./(fext0'*fext0+1000*eps));

    if true   % (15 Oct 2026) rBlocks always calculated
        Nrm=sqrt(fext0'*fext0+1000*eps);
        rBlocks=full([norm(frhs(1:N)) norm(frhs(N+1:2*N)) norm(frhs(2*N+1:3*N)) norm(grhs)]/Nrm);
        rBlocks=rBlocks .^2; 
    else
        rBlocks=[];
    end

end

end
