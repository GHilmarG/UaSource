function [b,h,GF]=Calc_b_From_sBS(CtrlVar,MUA,s,B,S,rho,rhow,GF,b0)


narginchk(7,9)
nargoutchk(3,3)


%% Calculates b from s, B, S and rho and rhow.
%
% GF and b0 are optional.
%
% GF and b0 are initial guesses for GF and b.
%
% Sets:
%
%   b=B over grounded areas.
%   b=(rho.*s-rhow.*S)./(rho-rhow)  over floating areas.
%
%
%  This is a bit dangerous function to use. This will not conserve thickness over the floating areas!
%
%
% This may make sense to do if one has measurments of the surface elevation and independent estimates of
% the grounding line, e.g. the floating mask.
%
% Because the floating mask depents on b through h, this is a non-linear
% problem.
%
% This is solved using the NR method. Usually only one single NR iteration is
% required.
%%

%% get a rough and a reasonable initial estimate for b

hf=rhow*(S-B)./rho ;
if nargin< 9  || isempty(b0)
    
    if nargin < 8 || isempty(GF)
        h=s-B ;  % for the purpose of calculating the floating mask, set b = B
        G = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);
    else
        G=GF.node;
    end
    
    b0 =  G.*B + (1-G).*(rho.*s-rhow.*S)./(rho-rhow) ;
end
%%

b=b0;
h=s-b;

% iteration
ItMax=30 ; tol=100*eps ;  Err=100*tol ; I=0 ;
ErrVector=zeros(ItMax,1)+NaN ;

while I< ItMax && Err > tol
    I=I+1;
    
    G = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);  % 1
    dGdb=-DiracDelta(CtrlVar.kH,h-hf,CtrlVar.Hh0) ;
    
    F0=    b - G.*B - (1-G).*(rho.*s-rhow.*S)./(rho-rhow) ;
    dFdb = 1 - dGdb.* (B -  (rho.*s-rhow.*S)./(rho-rhow)) ;
    
    db= -F0./dFdb ;
    
    b=b+db ;
    h=s-b ;
    
    F1 =    b - G.*B - (1-G).*(rho.*s-rhow.*S)./(rho-rhow) ;
    ErrLast=Err;
    Err=sum(F1.^2)/2 ;
    
    %fprintf('\t %i : \t %g \t %g \t %g \n ',I,max(abs(db)),Err,Err/ErrLast)
    
    ErrVector(I)=Err ;
    
end

GF.node = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);

% figure ; semilogy([1:30],ErrVector,'-or')


if I==ItMax   % if the NR iteration above, taking a blind NR step does not work, just
    % hand this over the matlab opt.
    % Why not do so right away? Because the above options is based on
    % my experience always faster if it converges (fminunc is very reluctant to take
    % lare steps, and apparantly does not take a full NR step...?!)
    
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