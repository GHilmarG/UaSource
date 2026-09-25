






function [dbdB,dhdB,dGdB,dFdb,d2bdB2,d2hdB2]=dGeometrydB(CtrlVar,s,S,B,b,rho,rhow)

%% Derivatives of the ice geometry with respect to the bedrock, B.
%
% Returns
%
%   db/dB , dh/dB , dG/dB
%
% for the geometrical closure in which the upper ice surface, s, the ocean surface, S, and the densities are held fixed,
% and b and h are calculated from s, S and B. This is the closure solved by
%
%   Calc_bh_From_sBS.m
%
% The closure is
%
% $$ b = \mathcal{G} B + (1-\mathcal{G}) b_s , \qquad h=s-b , \qquad \mathcal{G}=\mathcal{H}(h-h_f) $$
%
% with
%
% $$ b_s = \frac{\rho s - \rho_o S}{\rho-\rho_o} $$
%
% the base of freely floating ice, and
%
% $$ h_f = \frac{\rho_o (S-B)}{\rho} $$
%
% the flotation thickness.
%
%% Why this is an implicit problem
%
% The closure is implicit in b, because the floating mask depends on h=s-b. Writing
%
% $$ F_0(b,B) = b - \mathcal{G}(h-h_f) \, B - (1-\mathcal{G}(h-h_f)) \, b_s = 0 $$
%
% and using
%
% $$ \partial (h-h_f)/\partial b = -1 , \qquad \partial (h-h_f) / \partial B = \rho_o/\rho $$
%
% we get
%
% $$ \partial F_0/\partial b = 1 + \delta \, (B-b_s) $$
%
% $$ \partial F_0/\partial B = -\delta \, (\rho_o/\rho) \,(B-b_s) - \mathcal{G} $$
%
% where $\delta$ is the derivative of the (smoothed) Heaviside function. Note that the first of these is exactly the
% Newton Jacobian, dFdb, used within Calc_bh_From_sBS.m
%
% The implicit function theorem then gives
%
% $$ \frac{\partial b}{\partial B} = \frac{\mathcal{G}+ (\rho_o/\rho) \, \Lambda}{1+\Lambda} , \qquad \Lambda := \delta \,(B-b_s) $$
%
% $$ \frac{\partial h}{\partial B} = -\frac{\partial b}{\partial B} , \qquad \frac{\partial s}{\partial B}=0 $$
%
% and, using the total derivative of the mask,
%
% $$ \frac{\partial \mathcal{G}}{\partial B} = \delta \left ( \frac{\rho_o}{\rho} - \frac{\partial b}{\partial B} \right ) $$
%
%% Second derivative
%
% Differentiating the expression for db/dB once more with respect to B, and writing
%
% $$ \kappa = \rho_o/\rho , \qquad W = B-b_s , \qquad \Lambda = \delta W , \qquad \mu = \kappa - b' = \frac{d \Delta h}{dB} $$
%
% and using
%
% $$ \frac{d \mathcal{G}}{dB} = \delta \mu , \qquad \frac{d \Lambda}{dB} = \delta' \mu W + \delta , \qquad \mathcal{G}+\kappa \Lambda = b' (1+\Lambda) $$
%
% the algebra collapses to
%
% $$ b'' = \frac{ 2 \delta \mu + \delta' \, W \, \mu^2 }{1+\Lambda} , \qquad h'' = -b'' $$
%
% where $\delta'$ is the derivative of the smoothed delta with respect to its argument. For the tanh form used here
%
% $$ \delta = 2 k \, \mathcal{G} (1-\mathcal{G}) \qquad \Rightarrow \qquad \delta' = 2 k \, \delta \, (1-2\mathcal{G}) $$
%
% so no further special function is required.
%
% The second derivative is negligible away from the grounding line but is of order unity within the transition zone,
% where it is comparable in magnitude to db/dB itself. It is needed for the Hessian, where it multiplies the first
% derivative of the cost function with respect to h, and it contributes a diagonal term because the closure is applied
% independently at each node.
%
%% Important
%
% db/dB is NOT equal to the floating mask. Away from the grounding line it reduces to it (1 where grounded, 0 where
% afloat), but within the transition zone it exceeds unity on the grounded side (tending to rho_o/rho), and becomes
% NEGATIVE on the floating side. Both are genuine grounding-line migration effects: increasing B lowers h_f, which
% grounds the ice further, which pulls b from b_s towards B.
%
%% Where to evaluate
%
% The closure is applied at the nodes within Calc_bh_From_sBS.m, so b and h are ordinary finite-element fields whose
% nodal values depend on the nodal values of B. Therefore
%
% $$ \frac{\partial b(x)}{\partial B_j} = \left ( \frac{\partial b}{\partial B} \right )_j \, \phi_j(x) $$
%
% with the derivative evaluated AT NODE j. The quantities returned here are therefore nodal coefficients, and must not
% be interpolated onto the integration points. (They can be evaluated on any array, but if they are to be used as
% derivatives with respect to the nodal values of B, the inputs must be nodal.)
%
%% Inputs
%
% b must be the converged solution of the closure for the given B, i.e. the b returned by Calc_bh_From_sBS.m. The
% expressions above are exact only where F_0=0.
%
%  see also: Calc_bh_From_sBS.m, dIdBq.m
%
%%

narginchk(7,7)
nargoutchk(1,6)

bs = (rho.*s-rhow.*S)./(rho-rhow) ;    % base of freely floating ice
hf = rhow.*(S-B)./rho ;                % flotation thickness
Dh = (s-b)-hf ;                        % h-h_f

G     = HeavisideApprox(CtrlVar.kH,Dh,CtrlVar.Hh0) ;
Delta = DiracDelta(CtrlVar.kH,Dh,CtrlVar.Hh0) ;

Lambda = Delta.*(B-bs) ;
dFdb   = 1+Lambda ;                    % identical to dFdb within Calc_bh_From_sBS.m

if any(dFdb<=0)
    warning("dGeometrydB:SingularClosure",...
        "The geometrical closure b=b(B) is singular, or close to being so. min(dF0/db)=%g \n",min(dFdb))
end

dbdB = (G+(rhow./rho).*Lambda)./dFdb ;

if nargout>1
    dhdB = -dbdB ;
end

if nargout>2
    dGdB = Delta.*(rhow./rho-dbdB) ;
end

if nargout>4

    mu     = rhow./rho - dbdB ;                          % d(Dh)/dB
    dDelta = 2*CtrlVar.kH.*Delta.*(1-2*G) ;              % derivative of the smoothed delta

    d2bdB2 = ( 2*Delta.*mu + dDelta.*(B-bs).*mu.^2 )./dFdb ;

end

if nargout>5
    d2hdB2 = -d2bdB2 ;
end

%% Optional finite-difference test
% Set CtrlVar.TestGeometryDerivatives=true to compare db/dB against a central finite difference taken through
% Calc_bh_From_sBS.m itself.

if isfield(CtrlVar,"TestGeometryDerivatives") && isequal(CtrlVar.TestGeometryDerivatives,true)
    FiniteDifferenceTest(CtrlVar,s,S,B,rho,rhow,dbdB)
end

end




function FiniteDifferenceTest(CtrlVar,s,S,B,rho,rhow,dbdB)

CtrlVarTest=CtrlVar;

% Calc_bh_From_sBS only uses MUA within a diagnostic plotting branch, which is switched off here.
if ~isfield(CtrlVarTest,"MapOldToNew") || ~isfield(CtrlVarTest.MapOldToNew,"Test")
    CtrlVarTest.MapOldToNew.Test=false;
end

% The mask varies over a scale of about 1/kH, so the step must be small compared with that.
dB=1e-3/max(CtrlVar.kH,1e-10) ;

bPlus =Calc_bh_From_sBS(CtrlVarTest,[],s,B+dB,S,rho,rhow);
bMinus=Calc_bh_From_sBS(CtrlVarTest,[],s,B-dB,S,rho,rhow);

dbdB_FD=(bPlus-bMinus)/(2*dB);

Diff=norm(dbdB-dbdB_FD)/(norm(dbdB_FD)+eps);

fprintf("dGeometrydB: normalized norm of difference between analytical and FD db/dB is %g \n",Diff)

end
