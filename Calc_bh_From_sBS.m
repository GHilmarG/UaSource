function [b,h,GF]=Calc_bh_From_sBS(CtrlVar,MUA,s,B,S,rho,rhow)


narginchk(7,7)
nargoutchk(1,3)


%% 
%
% *Calculates b and h from s, B, S and rho and rhow.*
%
% GF and b0 are optional.
%
% GF and b0 are initial guesses for GF and b.
%
% Sets 
%
% $b=B$
%
% over grounded areas, and 
%
% $b=(\rho s-\rho_w S)/(\rho-\rho_w)$  
%
% over floating areas. 
%
% On return $s=b+h$.
%
% Note: This will not conserve thickness.
%
% Because the floating mask depents on b through h, this is a non-linear
% problem.
%
% Solved using the NR method. Usually only one single NR iteration is
% required.
%
%
% GF and b0   : (optional) initial guess for GF and b0
%
% MUA         : also optional and not currently used.
%
% Example:
%
%       b=Calc_bh_From_sBS(CtrlVar,[],s,B,S,rho,rhow)
%
%
%%

% get a rough and a reasonable initial estimate for b if none is provided
% The lower surface b is 
%
%
%   b=max( B , (rhow S - rho s)/(rhow-rho) ) 
%   where
%
%  h_f = rhow (S-B) / rho
%
%  b=s-h_f = 
hf=rhow*(S-B)./rho ;


b0 =  max(B,(rho.*s-rhow.*S)./(rho-rhow)) ; 

b=b0;
h=s-b;

% iteration
ItMax=30 ; tol=100*eps ;  J=Inf ; I=0 ;
JVector=zeros(ItMax,1)+NaN ;

while I< ItMax && J > tol
    I=I+1;
    
    G = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);  % 1
    dGdb=-DiracDelta(CtrlVar.kH,h-hf,CtrlVar.Hh0) ;
    
    F0=    b - G.*B - (1-G).*(rho.*s-rhow.*S)./(rho-rhow) ;
    dFdb = 1 - dGdb.* (B -  (rho.*s-rhow.*S)./(rho-rhow)) ;
    
    db= -F0./dFdb ;
    
    b=b+db ;
    h=s-b ;
    
    F1 =    b - G.*B - (1-G).*(rho.*s-rhow.*S)./(rho-rhow) ;
    
    JLast=J ;
    J=sum(F1.^2)/2 ;
    if CtrlVar.MapOldToNew.Test
        fprintf('\t %i : \t %g \t %g \t %g \n ',I,max(abs(db)),J,J/JLast)
    end
    
    JVector(I)=J ;
    
end

GF.node = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);

if CtrlVar.MapOldToNew.Test
    FindOrCreateFigure("Testing Calc_bh_From_sBs")  ; 
    semilogy([1:30],JVector,'-or')
end



if I==ItMax   % if the NR iteration above, taking a blind NR step does not work, just
    % hand this over the matlab opt.
    % Why not do so right away? Because the above options is based on
    % my experience always faster if it converges (fminunc is very reluctant to take
    % large steps, and apparantly does not take a full NR step...?!)
    
    options = optimoptions('fminunc','Algorithm','trust-region',...
        'SpecifyObjectiveGradient',true,'HessianFcn','objective',...
        'SubproblemAlgorithm','factorization','StepTolerance',1e-10,...
        'Display','iter');
    
    
    func=@(b) bFunc(b,CtrlVar,s,B,S,rho,rhow) ;
    b  = fminunc(func,b0,options) ;
    h=s-b;
    GF.node = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);
end
%%






end