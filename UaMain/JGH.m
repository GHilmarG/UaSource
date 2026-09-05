




function [J,dJdp,Hessian,JGHouts,F]=JGH(p,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint)


%%
%
% JGH: Returns the cost function (J), the gradient of the cost function with respect to p (dJdp), and the Hessian (ddJddp).
%
% The Hessian of the regularization term (R) can usually be calculated exactly, while the Hessian of the misfit/likelihood
% term (I), can not. However, one can come up with a educated guess for the Hessian of I with respect to C.
%
%
% Calculates objective function (J), gradient (dJdp, accurate), Hessian (guessed).
%
%
%%

persistent ubP vbP JGH1 JGH2 JGH3

narginchk(11,11)


%% some counters for how often JGH is called and with what number of arguments
if isempty(JGH1)
    JGH1=0;  % counter for 1-argument output, just cost function evaluation
    JGH2=0;  % counter for 2-argument output, cost and gradient 
    JGH3=0;  % counter for 3-argument output, cost, gradient and Hessian
end

if nargout==3

    JGH1=JGH1+1;
    JGH2=JGH2+1;
    JGH3=JGH3+1;

elseif nargout==2

    JGH1=JGH1+1;
    JGH2=JGH2+1;

elseif nargout==1

    JGH1=JGH1+1;

end

%% Create some internal flags indicating if cost, gradient or Hessian need to be calculated and returned. 

CtrlVar.Inverse.CalcGrad=false;
CtrlVar.Inverse.CalcGradI=false;
CtrlVar.Inverse.CalcGradR=false;


CtrlVar.Inverse.CalcHess=false;
CtrlVar.Inverse.CalcHessI=false;
CtrlVar.Inverse.CalcHessR=false;

if nargout==1
    dJdp=[] ; Hessian=[] ; JGHouts=[] ;
end

if nargout>=2  % always calculates the gradient if number of output arguments is 2 or larger
    CtrlVar.Inverse.CalcGrad=true;
    CtrlVar.Inverse.CalcGradI=true;
    CtrlVar.Inverse.CalcGradR=true;
end

% Calculate the Hessian provided:
%     1) CtrlVar.JGH.CalcHessian=true; 
% and 2) number of output arguments is larger or equal to 3
if nargout>= 3  &&  CtrlVar.JGH.CalcHessian 
    CtrlVar.Inverse.CalcHess=true;
    CtrlVar.Inverse.CalcHessI=true;
    CtrlVar.Inverse.CalcHessR=true;
else
    CtrlVar.Inverse.CalcHess=false;
    CtrlVar.Inverse.CalcHessI=false;
    CtrlVar.Inverse.CalcHessR=false;
    Hessian=[]; 
end

%% The function requires a solution of the forward model, often with just a slightly different input parameters. 
% Therefore, save previous velocity solution as persistent variables and use as initial staring point for next uv-solve.

if CtrlVar.Inverse.ResetPersistentVariables
    ubP=[];
    vbP=[];
end

if ~isempty(ubP)
    F.ub=ubP;
    F.vb=vbP;
end


%% Reflection?
% 
% if contains(CtrlVar.Inverse.MinimisationMethod,"Ua")
%     % pub and plb are enforced by the MATLAB optimization toolbox in a different way
%     if CtrlVar.ReflectiveTransformation
% 
%         if ~isempty(pub)  && isempty(plb)
% 
%             iu=p>pub;
%             p(iu)=pub(iu)+2*p(iu) ; % so if we had p(il)=plb(il) we get p(il)=plb(il)-2*p(il)=plb(il)
% 
%         elseif isempty(pub)  && ~isempty(plb)
% 
%             il=p>plb;
%             p(il)=pub(il)+2*p(il) ; % so if we had p(il)=plb(il) we get p(il)=plb(il)-2*p(il)=plb(il)
% 
%         elseif ~isempty(pub)  && ~isempty(plb)
% 
%             %%
%             % pub=[10 8]; plb=[1 2] ; p=[1 9] ;
% 
%             d=pub-plb;
%             t=mod(p-plb,2*d);
%             p=plb+min(t,2*d-t);
%             %%
% 
%         end
%     else
% 
%         p=kk_proj(p,pub,plb);
%     end
% 
% end
% 
%%




% The vector p contains the variables for which the inversion is being performed. So if the inversion is done over log(c)
% only, then p=log(C). And if the inversion is done over A, B and C then p=[A;B;C].

% Populate F with the current values in the vector p ahead of a call the the forward model.
F=p2F(CtrlVar,MUA,p,F,Meas,Priors);  % This maps from the vector p to the field variables F
                                     % p is the vector of control variables. This is some combination of A, B and C 
                                     % When inverting for A, we have p=[log10(AGlen]
                                     % When inverting for C, we have p=[log10(C)]
                                     % When inverting for A and C we have p=[log10(AGlen);log10(C)]
                                     % and so on
if anynan(F.C)
    error("JGH:Cnan","nan in C")
end
if anynan(F.AGlen)
    error("JGH:Anan","nan in A")
end
if anynan(F.B)
    error("JGH:Bnan","nan in B")
end

if any(F.C<0)
    error("JGH:Cneg","negative C values")
end

if any(F.AGlen<0)
    error("JGH:Cneg","negative A values")
end


%% Forward model solution
[~,~,F,l,dFduv]= uv([],[],CtrlVar,MUA,BCs,F,l);

if contains(CtrlVar.Inverse.Measurements,"-dhdt-")
    % If dh/dt is included as measurements, I need that calculated dh/dt, which in turn requires the mass-balance, a, as well
    if isempty(F.as) || isempty(F.ab)
        [~,F]=GetMassBalance([],CtrlVar,MUA,F);
    end
    [~,F.dhdt]=dhdtExplicit([],CtrlVar,MUA,F,BCs) ;
end

%%
% The cost function, J), is split into a misfit (I) and a regularization term (R). These usually consist of further
% terms.
%
% Get the I and R terms, and the gradients if required.
if nargout==1
    R=Regularisation(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint) ;
    I=Misfit(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,dFduv) ;
else
    [R,dRdp,ddRddp,RegOuts]=Regularisation(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint) ;
    [I,dIdp,ddIddp,MisfitOuts]=Misfit(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,dFduv) ;
end


if  CtrlVar.Inverse.CalcGrad  % gradient needed
    dJdp=dRdp+dIdp;
end

if CtrlVar.Inverse.CalcHess  % Hessian needed
    if isempty(ddIddp)
        Hessian=ddRddp;
    else
        Hessian=ddRddp+ddIddp;
    end
end








if F.solution=="-uv-"
    % To speed up the forward solve, the previous solution is stored locally and then used as a starting value in next
    % calculation. The idea is that usually the parameter vector (p) only changes slightly form one inverse iteration to the
    % next, so the (u,v) solution is likely to be similar to the previously calculated one.
    ubP=F.ub;
    vbP=F.vb;
else
    warning('JGH:returnsNaN',' uv solution did not converge. Returning NaN in cost function.\n ') ;
    ubP=[];
    vbP=[];
    I=NaN;
    R=NaN ;
    dJdp=p*0+NaN;
    MisfitOuts.I=NaN;
end

J=R+I;

if J < 0
    fprintf("J less that zero!! \n")
end


if nargout>3  % additional information needed as output
    JGHouts.dRdp=dRdp;
    JGHouts.dIdp=dIdp;
    JGHouts.ddIdpp=ddIddp;
    JGHouts.ddRdpp=ddRddp;
    JGHouts.RegOuts=RegOuts;
    JGHouts.MisfitOuts=MisfitOuts;
else
    JGHouts=[];
end

if  CtrlVar.Inverse.MinimisationMethod=="-MatlabOptimization-HessianVectorProduct-"
    % This is when using trust-region-reflective and Hessian-Vector-Product
    Hessian=p ;

end

if isnan(J)
    warning("JGH:ObjectivFunctionIsNaN","objective function is nan")
end


end

