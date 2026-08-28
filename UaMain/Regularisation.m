


function [R,dRdp,ddRdpp,RegOuts]=Regularisation(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint)

narginchk(8,8)


%%
% Calculates the regularization term R, and the gradient and the Hessian of R with respect to p.
%
%
%
% This is a fairly simple thing to do as the regularization term is an explicit function of p, and the Hessian calculation can
% be done exactly.
%
% However, there are quite a few cases to consider...
%
% *Regularisation Terms*
%
% The regularization terms can be thought of as the Priors in Bayesian context, provided they result in a valid covariance.
%
% Here the regularization terms are (mostly) on the form:
%
%
% $$ R = \int_{\mathcal{A}} \left ( \frac{1}{2} ( g f^2 +  k_x (\partial_x  f)^2 +  k_y (\partial_y  f)^2 )  \right ) \, \mathrm{d}\mathcal{A}   $$
%
% where $f$ is the deviation of the model parameter $p$ from the prior value $\tilde{p}$, i.e.
%
% $$f=p-\tilde{p}$$
%
% Taking the variation of $I$ results in the FE-equations
%
% $$ \delta R =  \langle g f | \delta f \rangle  + \langle k_x \partial_x f | \partial_x \delta f \rangle +  \langle k_y \partial_y f_y|\partial_y \delta f \rangle $$
%
% If we set $g=\gamma_a $ and $k_x=k_y = \gamma_s$ where $\gamma_a$ and $\gamma_s$ are constants, the FE discretized system is
%
%
% $$\delta R =  ( \gamma_a  \mathbf{M} + \gamma_s (\mathbf{D}_x + \mathbf{D}_y ) ) \, \mathbf{f} $$
%
% and
%
% $$R =  \frac{1}{2} \mathbf{f}^T  ( \gamma_a  \mathbf{M} + \gamma_s (\mathbf{D}_x +\mathbf{D}_y ) ) \mathbf{f} $$
%
% The discretized precision matrix $\mathbf{P}$ is therefore
%
% $$ \mathbf{P}= \frac{1}{2} ( \gamma_a  \mathbf{M} + \gamma_s (\mathbf{D}_x +\mathbf{D}_y ) ) $$
%
% The precision matrix relates to our estimated covariance between $p$ and $\tilde{p}$. We may have some additional
% (direct) information about $p$ in the form of some further measurements of $p$ not included in our prior estimate
% $\tilde{p}$. If our prior for $p$ is $P(p)$ before we obtain the data $\hat{p}$, and $P(\hat{p} | p)$ is the likelihood of
% the data given $p$, then updated estimate for $p$ after having seen the data is
%
% $$ P(p|\hat{p}) \propto P(\hat{p}|p) P(p) $$
%
%
%
% The precision matrix above can be thought of as the inverse covariance of the noise free latent (i.e. unobserved) variable.
%
%
%
% For each variable (A,B,C) the regularization term typically has the form
%
%   R= ( ga.^2* Ra + gs.^2 * Rs ) /(2*Area)
%
% where ga and gs are regularization parameters
%
% and
%
%  Ra = dp'*M*dp
%
%  Rs= dp'*(Dxx+ Dyy)*dp
%
% where dp = p-pPrior, and M is the mass matrix and Dxx and Dyy the stiffness matrices.
%
% *Misfit Terms*
%
%
% We might also have direct measurements of the latent fields (A, B, C). For example, when we have direct measurements of B.
% This are really misfit terms, but I include those here because they are also explicit functions of the latent variables.
%
% B 'misfit' term has the form:
%
% $$ \int (B-B_{\mathrm{meas}}) \, \mathcal{E}^{-1} \, (B-B_{\mathrm{meas}}) \, dA $$
%
% where
%
% $\mathcal{E}$
%
% is the error covariance.
%
% By defining:
%
% $$ \tilde{B} = (B-B_{\mathrm{meas}})/B_{\epsilon} $$
%
% This simply has the form
%
% $$ \tilde{B} \, M \tilde{B} $$
%
% *Terms involving functions of the latent variables*
%
% If we have a term involving a function of a variable, e.g.
%
% $$f(B)$$
%
% with
%
% $$
% R =\frac{1}{2} \int (f(B))^2 \, \mathrm{d}x
% = \frac{1}{2} \int (f(B)_i \phi_i(x) ) \, ( f(B)_j \phi_j(x) ) \, dx
% = \frac{1}{2} f_i(B) \left  (  \int \phi_i(x) ) \, \phi_j(x) \, dx \right ) \, f_j(B)
% = \frac{1}{2} f_i(B) \, M_{ij} \, f_j(B)
% $$
%
% where, for example, we could have
%
% $$ f(B)=\mathrm{SoftPlus}(B-s+h_{\mathrm{min}}) $$
%
%
% $$ f(B) = \sum_{k=1}^n \, f_k(B_k) \, \phi_k(x) $$
%
% Denote the derivative with $g$ , i.e.
%
% $$ g(B) := \partial f(B)/\partial B = \sum_{i=1}^n \, \frac{\partial f_k}{\partial B} \, \phi_k(x) = \sum_{k=1}^n \, g_k \, \phi_k(x) $$
%
% Note that in the case where, at every node
%
% $$g_k=\frac{\partial f_k}{\partial B}=1 $$
%
% for every $g_k$, we have
%
% $$
% g(B) := \partial f(B)/\partial B
% = \sum_{k=1}^n \, \frac{\partial f_k}{\partial B} \, \phi_k(x)
% = \sum_{k=1}^n \, g_k \, \phi_k(x)
% = \sum_{k=1}^n \, 1 \, \phi_k(x)  = 1
% $$
%
% due to the partition of unity property of the finite element form functions
%
% $$\sum_{k=1}^n \, \phi_k(x)  = 1 $$
%
%
% In general
%
% $$
% \delta_B R(B) = \lim_{\epsilon \to 0} \, \frac{1}{2} \frac{d}{d \epsilon} \int  (f(B(x)+\epsilon \phi_q(x)))^2 \, dx
% =\int f(B) \, g(B)\, \phi_q  \; dx
% = \int \; f_i \phi_i(x) \; g_k \phi_k(x) \; \phi_q(x)  \; dx
% $$
%
% In the above expression, we have taken the derivatives of $f(B)$ at the nodes and then interpolated those nodal values to the
% integration points. We could alternative interpolate $B$ to the integration points, and then take the derivative of $f(B)$
% with respect of $B$ at the location of the Integration points.
%
% Doing so results in:
%
% $$
% \delta_B R(B)
% = \int \; f_i \phi_i(x) \; \frac{\partial f(B_k \phi_k(x))}{\partial B}  \; \phi_q(x)  \; dx
% $$
%
% The term
%
% $$ \frac{\partial f(B_k \phi_k(x))}{\partial B}$$
%
% is a scalar, i.e. this derivative at a given location $x$
%
% Generally, when taking several derivatives and making use of the chain rule, it is best to evaluation derivative and other
% functions of the primary variables at the integration points.
%
% The Hessian is then
%
% $$
% H = \delta_{BB} I(B)
% = \int \; \frac{\partial f(B_i \phi_i(x)) }{\partial B} \, \phi_r(x) \; \frac{\partial f(B_k \phi_k(x))}{\partial B}  \; \phi_q(x)  \; dx
% + \int \; f_i \phi_i(x)  \; \frac{\partial^2 f(B_k \phi_k(x))}{\partial B \, \partial B}  \; \phi_r(x) \, \phi_q(x)  \; dx
% $$
%
% Note that for
%
% $$ f(B)= B = B_i \phi_i(x) $$
%
% in which case
%
% $$ f_i=B_I $$
%
% and
%
% $$f^{\prime} =1 $$
%
% and
%
% $$f^{\prime\,\prime} =0  $$
%
%
% we have
%
% $$ \frac{\partial p}{\partial B_j} = \delta_{ij} $$
%
%
% and we arrive at
%
% $$
% \delta_B R(B) = \lim_{\epsilon \to 0} \, \frac{1}{2} \frac{d}{d \epsilon} \int  (f(B(x)+\epsilon \phi_q(x)))^2 \, dx
% = f_i \; \left ( \int \phi_i(x) \, \phi_j(x) \, dx \right ) \frac{\partial f(B)}{\partial B_j}
% = f_i M_{ij} \phi_j \, \delta_{ij} = B_i M_{ij} = M_{ji} B_j
% $$
%
% or we then have (as listed above)
%
% $$ R= \frac{1}{2} \mathbf{B}^T \mathbf{M} \mathbf{B} $$
%
% $$\delta_B R(B) = \mathbf{M} \mathbf{B} $$
%
% and
%
% $$ H= \mathbf{M} $$
%
%
%%


if nargout > 3
    RegOuts=[];
    RegOuts.RAs=nan  ; RegOuts.RAa=nan;
    RegOuts.RCs=nan  ; RegOuts.RCa=nan;

end
%%

% Add up field by field

RA=0;
dRdA=[];
ddRdAA=[];

RC=0;
dRdC=[];
ddRdCC=[];


Area=MUA.Area;


%%
% I start by defining dpX, gsX and gaX, where X is either C or A and
%  dpX=X-X_{Prior}
%  gsX and gaX are the slope and amplitude regularization pre-factors
%

% C
if contains(lower(CtrlVar.Inverse.InvertFor),'c')  % this includes both c and logc inversion

    isC=1;
    if contains(lower(CtrlVar.Inverse.Regularize.Field),'logc')

        % regularize log10(C)
        dpC=log10(F.C)-log10(Priors.C);
        pPriorCovC=Priors.CovC;
        gsC=CtrlVar.Inverse.Regularize.logC.gs;
        gaC=CtrlVar.Inverse.Regularize.logC.ga;

        if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
            dCfactor=1;
        else
            dCfactor=1./F.C/log(10); % d (dpC)/dC= 1/(log(10) C)
        end

    else
        % regularize C
        dpC=F.C-Priors.C;
        pPriorCovC=Priors.CovC;
        gsC=CtrlVar.Inverse.Regularize.C.gs;
        gaC=CtrlVar.Inverse.Regularize.C.ga;
        if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
            dCfactor=log(10)*F.C;   % gradient must be with respect to logC, but regularisation is on C
        else
            dCfactor=1;
        end
    end
else

    isC=0;
    dpC=0;
    dCfactor=0;
    pPriorCovC=1;
    gsC=0;
    gaC=0;

end

% AGlen
if contains(lower(CtrlVar.Inverse.InvertFor),'aglen')

    isA=1;
    if contains(lower(CtrlVar.Inverse.Regularize.Field),'logaglen')

        % regularize log10(AGlen)

        dpA=log10(F.AGlen)-log10(Priors.AGlen);
        pPriorCovA=Priors.CovAGlen;
        gsA=CtrlVar.Inverse.Regularize.logAGlen.gs;
        gaA=CtrlVar.Inverse.Regularize.logAGlen.ga;
        if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            dAfactor=1;
        else
            dAfactor=1./F.AGlen/log(10);   % gradient must be with respect to A, but regularization is on logA
        end

    else % regularize A

        dpA=F.AGlen-Priors.AGlen;
        pPriorCovA=Priors.CovAGlen;
        gsA=CtrlVar.Inverse.Regularize.AGlen.gs;
        gaA=CtrlVar.Inverse.Regularize.AGlen.ga;
        if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            dAfactor=log(10)*F.AGlen;   % gradient must be with respect to logA, but regularization is on A
        else
            dAfactor=1;
        end
    end
else

    isA=0;
    dpA=0;
    dAfactor=0;
    pPriorCovA=1;
    gsA=0;
    gaA=0;

end

% b
if contains(CtrlVar.Inverse.InvertFor,'-b-')
    error('fdsa')
    isb=1;
    dpb=F.b-Priors.b;
    pPriorCovb=Priors.Covb;
    gsb=CtrlVar.Inverse.Regularize.b.gs;
    gab=CtrlVar.Inverse.Regularize.b.ga;
    dbfactor=1;

else

    isb=0;
    dpb=0;
    dbfactor=0;
    pPriorCovb=1;
    gsb=0;
    gab=0;

end


% B
if contains(CtrlVar.Inverse.InvertFor,'-B-')

    isB=1;
    dpB=F.B-Priors.B;
    pPriorCovB=Priors.CovB;
    gsB=CtrlVar.Inverse.Regularize.B.gs;
    gaB=CtrlVar.Inverse.Regularize.B.ga;
    dBfactor=1;

else

    isB=0;
    dpB=0;
    dBfactor=0;
    pPriorCovB=1;
    gsB=0;
    gaB=0;

end

if ~isfield(MUA,'M') || isempty(MUA.M)
    MUA.M=MassMatrix2D1dof(MUA);
end

if ~isfield(MUA,'Dxx') || isempty(MUA.Dxx)
    [MUA.Dxx,MUA.Dyy]=StiffnessMatrix2D1dof(MUA);
end

M=MUA.M;
Dxx=MUA.Dxx;
Dyy=MUA.Dyy;



%% Now dpX, gsX and gaX have all be defined

%% Now defining R, dRdp, ddRddp


% the expression for the prior, is
%
% $$-\log P(B) = \frac{1}{2}(B-B_{prior})^{T} Q (B-B_{prior}) + \frac{1}{2}\log\left|Q^{-1}\right| + \text{const} $$

if isA



    alphaMatern=CtrlVar.Inverse.Matern.logAGlen.alpha;
    kappaMatern=CtrlVar.Inverse.Matern.logAGlen.kappa;
    tauMatern=CtrlVar.Inverse.Matern.logAGlen.tau;
    gsA=CtrlVar.Inverse.Regularize.logAGlen.gs;
    gaA=CtrlVar.Inverse.Regularize.logAGlen.ga;

    QA=PrecisionMatrixMatern(MUA,alphaMatern,kappaMatern,tauMatern,gaA,gsA,CtrlVar.Inverse.Methodology);
    RA=0.5*dpA'*QA*dpA;           % costs function term
    dRdA=(QA*dpA).*dAfactor;      % derivative, accounting for possible log

    if nargout > 3  % this is here for a possible info, and useful when calculating L curves
  
        RegOuts.RAs=nan  ; RegOuts.RAa=nan;  % this is no longer working using Matern, I need to re-introduce this somehow, to be done
    end



    if  contains(CtrlVar.Inverse.MinimisationMethod,"HessianFiniteDifferences")
        N=MUA.Nnodes;
        ddRdAA=sparse(N,N);
    elseif contains(CtrlVar.Inverse.MinimisationMethod,"Hessian")

        if contains(CtrlVar.Inverse.Hessian,"RHA=E")
            ddRdAA=QA.*dAfactor;
        elseif contains(CtrlVar.Inverse.Hessian,"RHA=M")
            ddRdAA=MUA.M/MUA.Area;
        elseif contains(CtrlVar.Inverse.Hessian,"RHA=I") || contains(CtrlVar.Inverse.Hessian,"RHA=1")
            N=MUA.Nnodes;
            ddRdAA=speye(N,N);
        elseif contains(CtrlVar.Inverse.Hessian,"RHA=0") || contains(CtrlVar.Inverse.Hessian,"RHA=O")
            N=MUA.Nnodes;
            ddRdAA=sparse(N,N);
        else
            fprintf(" CtrlVar.Inverse.Hessian=%s, is incorrect.\n",CtrlVar.Inverse.Hessian)
            error("Regularisation:IncorrectInputs"," case not found ")
        end
    end

end

if isC



    %QC=(gsC.^2.*(Dxx+Dyy)+gaC.^2.*M)/Area;

    alphaMatern=CtrlVar.Inverse.Matern.logC.alpha;
    kappaMatern=CtrlVar.Inverse.Matern.logC.kappa;
    tauMatern=CtrlVar.Inverse.Matern.logC.tau;
    gsC=CtrlVar.Inverse.Regularize.logC.gs;
    gaC=CtrlVar.Inverse.Regularize.logC.ga;

    QC=PrecisionMatrixMatern(MUA,alphaMatern,kappaMatern,tauMatern,gaC,gsC,CtrlVar.Inverse.Methodology);
    RC=0.5*dpC'*QC*dpC;           % costs function term
    dRdC=(QC*dpC).*dCfactor;      % derivative, accounting for possible log

    if nargout > 3  % this is here for a possible info, and useful when calculating L curves
    
        RegOuts.RCs=nan  ; RegOuts.RCa=nan; % this is no longer working using Matern, I need to re-introduce this somehow, to be done
    end

    if  contains(CtrlVar.Inverse.MinimisationMethod,"HessianFiniteDifferences")
        N=MUA.Nnodes;
        ddRdCC=sparse(N,N);
    elseif contains(CtrlVar.Inverse.MinimisationMethod,"Hessian")
        if contains(CtrlVar.Inverse.Hessian,"RHC=E")

            ddRdCC=QC.*dCfactor;
        elseif contains(CtrlVar.Inverse.Hessian,"RHC=M")
            ddRdCC=MUA.M/MUA.Area;
        elseif contains(CtrlVar.Inverse.Hessian,"RHC=I") || contains(CtrlVar.Inverse.Hessian,"RHC=1")
            N=MUA.Nnodes;
            ddRdCC=speye(N,N);
        elseif contains(CtrlVar.Inverse.Hessian,"RHC=O") || contains(CtrlVar.Inverse.Hessian,"RHC=0")
            N=MUA.Nnodes;
            ddRdCC=sparse(N,N);
        else

            fprintf(" CtrlVar.Inverse.Hessian=%s, is incorrect.\n",CtrlVar.Inverse.Hessian)
            error("Regularisation:IncorrectInputs"," case not found ")

        end
    end

end


if isb   %  b

    Nb=(gsb.^2.*(Dxx+Dyy)+gab.^2.*M)/Area;
    Rb=dpb'*Nb*dpb/2;
    dRdb=(Nb*dpb).*dbfactor;
    error('fdsa')
else
    Rb=0;
    dRdb=[];
end



if isB   %  B

    % The covariance of the prior consists of two terms:
    %
    % 1) The 'usual' large-scale correlation which here is a Marten covariance,
    % 2) A 'nugget' effect which is related to uncorrelated errors in the (direct) measurements of B.
    %
    %QA=0.5*(gsA.^2.*(Dxx+Dyy)+gaA.^2.*M)/Area; % This is the precision matrix

    alphaMatern=CtrlVar.Inverse.Matern.B.alpha;
    kappaMatern=CtrlVar.Inverse.Matern.B.kappa;
    tauMatern=CtrlVar.Inverse.Matern.B.tau;
    QB=PrecisionMatrixMatern(MUA,alphaMatern,kappaMatern,tauMatern,gaB,gsB,CtrlVar.Inverse.Methodology);
    %QB=(gsB.^2.*(Dxx+Dyy)+gaB.^2.*M)/Area;
    RB=dpB'*QB*dpB/2;               %       R: Regularisation term for B (a scalar)
    dRdB=(QB*dpB).*dBfactor;        %   dR/dB:  (a vector)
    ddRdBB=QB.*dBfactor;            % exact, or simply the correct, Hessian of the regularization term
    % To do: I could add "RHB=E" to CtrlVar.Inverse.Hessian. Right now I do the exact (E) Hessian evaluation here.


    if ~isempty(Meas.B)  &&  ~isempty(Meas.BCov)  &&    isdiag(Meas.BCov)

        % Adding a cost term giving the deviation of inverted B from direct measurements of B. This has the same form as a data
        % misfit term used for velocities and dh/dt. But here this is applied to the inverted field.
        %
        % It could be argued that this term should be added to the likelihood (i.e. the misfit term) but here this
        % distinction is simply rhetorical as these terms are all added up

        Berr=sqrt(spdiags(Meas.BCov));

        Bres=(F.B-Meas.B)./Berr;
        RBmeas=full(Bres'*MUA.M*Bres)/2/Area;
        dRdBmeas=(MUA.M*Bres)./Berr/Area;
        ddRdBmeasBmeas=(MUA.M)./Berr/Area;

        RB=RB+RBmeas;
        dRdB=dRdB+dRdBmeas;
        ddRdBB=ddRdBB+ddRdBmeasBmeas;

    end


    %%  Barrier term to push B solution away from min ice thickness, i.e. to discourage F.B being close to F.s/Meas.s#
    %
    % Idea:  Add a quadratic penalty in terms of min thickness violation.
    %
    % Thickness violation: F.s - F.B < hmin

    CtrlVar.Inverse.Penalty=false;

    if CtrlVar.Inverse.Penalty

        x=F.B - (F.s-20*CtrlVar.ThickMin);
        x0=zeros(MUA.Nnodes,1);
        k=0.1; a=5;  % 1/k is the softness and a the amplitude
        [Bbarr,dBarrdB,ddBbarrdBB]=JgHpenalty(CtrlVar,MUA,x,x0,k,a) ;

        Bbarr=Bbarr/Area;
        dBarrdB=dBarrdB/Area;
        ddBbarrdBB=ddBbarrdBB/Area;


        RB=RB+Bbarr;
        dRdB=dRdB+dBarrdB;
        ddRdBB=ddRdBB+ddBbarrdBB;

    end

    if  contains(CtrlVar.Inverse.MinimisationMethod,"HessianFiniteDifferences")
        N=MUA.Nnodes;
        ddRdBB=sparse(N,N);
    end
    %

    %%


else
    RB=0;
    dRdB=[];
    ddRdBB=[];
end


% if CtrlVar.Inverse.MinimisationMethod contains "Hessian", then the pre-multipler is simply I, so this has no effect.
dRdA=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,BCsAdjoint,CtrlVar.Inverse.AdjointGradient.UseBCs.A,dRdA);
dRdC=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,BCsAdjoint,CtrlVar.Inverse.AdjointGradient.UseBCs.C,dRdC);
dRdB=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,BCsAdjoint,CtrlVar.Inverse.AdjointGradient.UseBCs.B,dRdB);

R=RA+RB+RC;
dRdp=[dRdA;dRdB;dRdC];
ddRdpp=blkdiag(ddRdAA,ddRdBB,ddRdCC) ; % much faster



R=CtrlVar.Inverse.Regularize.Multiplier*R;
dRdp=CtrlVar.Inverse.Regularize.Multiplier*dRdp;
ddRdpp=CtrlVar.Inverse.Regularize.Multiplier*ddRdpp;

if nargout > 3
    RegOuts.R=R;
    RegOuts.dRdp=dRdp;
    RegOuts.ddRddp=ddRdpp;


    RegOuts.RAGlen=RA;
    RegOuts.dRdAGlen=dRdA;

    RegOuts.RC=RC;
    RegOuts.dRdC=dRdC;

    RegOuts.Rb=Rb;
    RegOuts.dRdb=dRdb;

    RegOuts.RB=RB;
    RegOuts.dRdB=dRdB;
end

if R< 0
    fprintf("Regularisation.m : R is negative \n")
end


end