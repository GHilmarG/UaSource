





function [b,h,GF]=Calc_bh_From_sBS(CtrlVar,MUA,s,B,S,rho,rhow,G0)


narginchk(7,8)
nargoutchk(1,3)

if nargin < 8 || isempty(G0)
    G0=nan;
end

if isstruct(G0)
    if isfield(G0,"node")
        G0=G0.node;
    end
end


%% Geometrical closure: calculates b and h from s, B, S, rho and rhow
%
%   [b,h,GF]=Calc_bh_From_sBS(CtrlVar,MUA,s,B,S,rho,rhow)
%   [b,h,GF]=Calc_bh_From_sBS(CtrlVar,MUA,s,B,S,rho,rhow,G0)
%
% The upper ice surface, $s$, the ocean surface, $S$, the bedrock, $B$, and the densities, $\rho$ and $\rho_o$ (rhow), are given, and are held fixed.
% The lower ice surface, $b$, and the ice thickness, $h$, are calculated. This is the geometrical closure used when $s$, $B$ and $S$ are the
% independent geometrical variables, as in a $B$ inversion (CtrlVar.Calculate.Geometry="bh-FROM-sBS"). It is the counterpart of
% Calc_bs_From_hBS.m, where $h$, $B$ and $S$ are the independent variables.
%
%% The closure
%
% Define the flotation thickness
%
% $$ h_f = \frac{\rho_o}{\rho} \, (S-B) , $$
%
% and the base of freely floating ice with upper surface $s$
%
% $$ b_s = \frac{\rho s - \rho_o S}{\rho-\rho_o} . $$
%
% The solution of the closure, $b_c$ and $h_c$, is then
%
% $$ b_c = \mathcal{G} \, B + (1-\mathcal{G}) \, b_s , \qquad h_c = s - b_c , \qquad \mathcal{G} = \mathcal{H}(h_c - h_f) , $$
%
% where $\mathcal{H}$ is the smoothed Heaviside function (HeavisideApprox.m). Hence $b_c=B$ where grounded ($\mathcal{G}=1$), and $b_c=b_s$ where
% floating ($\mathcal{G}=0$). Equivalently
%
% $$ h_c = \mathcal{G} \, (s-B) + (1-\mathcal{G}) \, (s-b_s) . $$
%
% Because the floating mask depends on $h_c$, and hence on $b_c$, this is a non-linear equation for $b_c$. It is solved using the Newton-Raphson
% method applied to
%
% $$ F_0(b) = b - \mathcal{G} \, B - (1-\mathcal{G}) \, b_s , \qquad \frac{\partial F_0}{\partial b} = 1 + \delta \, (B-b_s) , $$
%
% where $\delta = \mathcal{H}^{\prime}(h-h_f)$ is the derivative of the smoothed Heaviside function. Usually only one single iteration is required.
%
% Note: This does not conserve thickness.
%
%% Minimum thickness
%
% With $h_{\min}$=CtrlVar.ThickMin the thickness returned is not $h_c$, but
%
% $$ h = \Phi(h_c) = h_{\min} + w \, \ln \left( 1 + \exp \left( \frac{h_c-h_{\min}}{w} \right) \right) , \qquad b = s - h , \qquad \mathcal{G}_{\mathrm{out}} = \mathcal{H}(h-h_f) $$
%
% where the width, $w$, is CtrlVar.ThickMinWidthRelative times $h_{\min}$ (default 0.05). This is a smooth floor: $h>h_{\min}$ everywhere, and
% $h=h_c$ to machine precision where $h_c>h_{\min}+40 w$. At $h_c=h_{\min}$ the thickness is increased by $w \ln 2$ (0.035 m by default). The
% floor is defined, and its inverse and derivatives are calculated, in ThicknessFloor.m.
%
% The upper surface $s$ is held fixed, so $b+h=s$ also where the floor is active. At those nodes $b$ is therefore no longer the solution of the closure
% above. Because $h>h_{\min}$, geometry calculated here never triggers the reset of thicknesses below ThickMin in uv.m, which otherwise rebuilds $b$ and $s$
% at ALL nodes (using Calc_bs_From_hBS.m), a rebuild that is not consistent with the closure within the grounding-line transition zone.
%
% The derivatives of $b$ and $h$ with respect to $B$ are those of the closure, scaled by $\Phi^{\prime}$ (see dGeometrydB.m)
%
% $$ \frac{dh}{dB} = \Phi^{\prime}(h_c) \, \frac{dh_c}{dB} , \qquad \frac{db}{dB} = - \frac{dh}{dB} , \qquad \frac{d^2 h}{dB^2} = \Phi^{\prime\prime}(h_c) \left ( \frac{dh_c}{dB} \right )^2 + \Phi^{\prime}(h_c) \, \frac{d^2 h_c}{dB^2} $$
%
% The floor is switched off by setting CtrlVar.ThickMinWidthRelative=0, in which case the closure above is returned unchanged.
%
%% Relation to Calc_bs_From_hBS.m
%
% In a forward run the thickness, $h$, is the prognostic variable, $B$ is fixed, and $s$ evolves. The geometry is then calculated from $h$ and $B$ by
% Calc_bs_From_hBS.m, in which the floating base is calculated from the thickness,
%
% $$ b = \mathcal{G} \, B + (1-\mathcal{G}) \, ( S - \rho h/\rho_o ) , \qquad s = b + h , \qquad \mathcal{G} = \mathcal{H}(h-h_f) $$
%
% This leaves $h$ unchanged, and therefore conserves thickness.
%
% In a $B$ inversion, $s$ is data and $B$ is the parameter. A change in $B$ must then change $b$ and $h$ at fixed $s$, so there is no thickness to conserve, and the
% closure above, in which the floating base, $b_s$, is calculated from $s$, is the appropriate one.
%
% The two closures agree where $\mathcal{G}=0$ and $\mathcal{G}=1$, but they differ within the grounding-line transition zone, $0<\mathcal{G}<1$ (by about 1 m in
% thickness for CtrlVar.kH=1 when the geometry from one is passed through the other). This function is therefore not the exact inverse of Calc_bs_From_hBS.m within the
% transition zone. When geometry from a $B$ inversion is used in a forward run, $b$ and $s$ are recalculated from $h$ and $B$ using Calc_bs_From_hBS.m, and differ
% slightly, at the nodes within the transition zone, from those used in the inversion.
%
%% Inputs
%
% G0 can be either an initial guess for the nodal grounded/floating mask itself, e.g. GF.node or F.GF.node. But it can also be GF where GF.node is
% the grounded/floating mask.
%
% G0 is just an initial guess, and simply omitting it as an input is perfectly fine.
%
% MUA         : also optional and not currently used.
%
% Example:
%
%       b=Calc_bh_From_sBS(CtrlVar,[],s,B,S,rho,rhow)
%
%  see also: ThicknessFloor.m, dGeometrydB.m, Calc_bs_From_hBS.m, p2F.m
%
%%

% get a rough and a reasonable initial estimate for b
% The lower surface b is 
%
%
%   b=max( B , (rhow S - rho s)/(rhow-rho) ) 
%   where
%
%  h_f = rhow (S-B) / rho
%
%  b=s-h_f 






hf=rhow*(S-B)./rho ;

% b_f is b based on flotation condition. This is only a function of s and S and the densities, and as these do not change in
% the course of the iteration, b_f remains the same.
b_f=(rho.*s-rhow.*S)./(rho-rhow);  

% a rough initial estimate for b.
% For this initial iterate the Newton iteration will converge for \rho_o/\rho>1 
b0 =  max(B,b_f) ; 



b=b0;
h=s-b;

% iteration
ItMax=30 ; tol=1000*eps ;  J=Inf ; I=0 ;
JVector=zeros(ItMax,1)+NaN ;

while I < ItMax && J > tol
    I=I+1;

    if isnan(G0)
        G = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);  % 1
        dGdb=-DiracDelta(CtrlVar.kH,h-hf,CtrlVar.Hh0) ;
    else
        G=G0;
        dGdb=0;

    end

    F0=    b - G.*B - (1-G).*b_f ;
    dFdb = 1 - dGdb.* (B -  b_f) ;
    
    db= -F0./dFdb ;
    
    b=b+db ; % b is updated,
    h=s-b ;  % h is updated
    
    F1 =    b - G.*B - (1-G).*b_f ;
    
    JLast=J ;
%    J=sum(F1.^2)/2 ;
    J=norm(F1)^2/2/numel(F1) ;

    if CtrlVar.MapOldToNew.Test
        fprintf('\t %i : \t %g \t %g \t %g \n ',I,max(abs(db)),J,J/JLast)
    end
    
    JVector(I)=J ;

    if J< tol
        break
    end
    
end


% smooth floor on the thickness, h>=ThickMin, with s fixed (see the header)
[hFloor,~,~,~,FloorIsActive]=ThicknessFloor(CtrlVar,h,"forward") ;
if FloorIsActive
    h=hFloor ;
    b=s-h ;
end

GF.node = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);

if CtrlVar.MapOldToNew.Test
    FindOrCreateFigure("Testing Calc_bh_From_sBs")  ; 
    semilogy(1:30,JVector,'-or')
    xlabel("iterations",Interpreter="latex")
    ylabel("Cost function, $J$",Interpreter="latex")
    title("Calculating $b$ and $h$ from $s$, $S$, and $B$",Interpreter="latex")
    title(sprintf("Calculating $b$ and $h$ from $s$, $S$, and $B$ by minimizing \n $J=\\int (b-\\mathcal{G}B - (1-\\mathcal{G}) (\\rho s -\\rho_o S/(\\rho-\\rho_o))\\, \\mathrm{d}x \\, \\mathrm{d}y$\n with respect to $b$ "),Interpreter="latex")

    % f=gcf ; exportgraphics(f,'Calc_bh_from_sBS_Example.pdf')

end



if I==ItMax   % if the NR iteration above, taking a blind NR step does not work, just
    % hand this over the matlab opt.
    % Why not do so right away? Because the above options is based on
    % my experience always faster if it converges (fminunc is very reluctant to take
    % large steps, and apparently does not take a full NR step...?!)
    %
    % Also, the Newton iteration is guaranteed to converge from the first iterate above provided rho_o/rho > 1
    %
    
    warning("Calc_bh_From_SBS:NoConvergence","Calc_bh_from_sBS did not converge! \n")

    options = optimoptions('fminunc','Algorithm','trust-region',...
        'SpecifyObjectiveGradient',true,'HessianFcn','objective',...
        'SubproblemAlgorithm','factorization','StepTolerance',1e-10,...
        'Display','iter');
    
    
    func=@(b) bFunc(b,CtrlVar,s,B,S,rho,rhow) ;
    b  = fminunc(func,b0,options) ;
    h=s-b;
    [hFloor,~,~,~,FloorIsActive]=ThicknessFloor(CtrlVar,h,"forward") ;
    if FloorIsActive
        h=hFloor ;
        b=s-h ;
    end
    GF.node = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);
end
%%






end