function [b,h,GF]=Calc_bh_From_sBS(CtrlVar,MUA,s,B,S,rho,rhow)


narginchk(7,7)
nargoutchk(1,3)


%% 
%
% *Calculates b and h from s, B, S and rho and rhow.*
%
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
% Because the floating mask depends on b through h, this is a non-linear
% problem.
%
% Solved using the NR method. Usually only one single NR iteration is
% required.
%
%
%
% MUA         : also optional and not currently used.
%
% Example:
%
%       b=Calc_bh_From_sBS(CtrlVar,[],s,B,S,rho,rhow)
%
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


b0 =  max(B,(rho.*s-rhow.*S)./(rho-rhow)) ; % a rough initial estimate for b

b=b0;
h=s-b;

% iteration
ItMax=30 ; tol=1000*eps ;  J=Inf ; I=0 ;
JVector=zeros(ItMax,1)+NaN ;

while I < ItMax && J > tol
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

    if J< tol
        break
    end
    
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
    % large steps, and apparantly does not take a full NR step...?!)
    
    warning("Calc_bh_From_SBS:NoConvergence","Calc_bh_from_sBS did not converge! \n")

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