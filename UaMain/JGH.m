




function [J,dJdp,Hessian,JGHouts,F,RunInfo]=JGH(p,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)


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

narginchk(14,14)
CtrlVar.nargoutJGH=nargout;

if isempty(JGH1)
    JGH1=0;
    JGH2=0;
    JGH3=0;
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

fprintf("JGH: 1 %i \t 2 %i \t 3 %i \n",JGH1,JGH2,JGH3)



if nargout==1
    CtrlVar.Inverse.CalcGradI=false;
    CtrlVar.Inverse.CalcGradR=false;
    dJdp=[] ; Hessian=[] ; JGHouts=[] ;
else
    CtrlVar.Inverse.CalcGradI=true;
    CtrlVar.Inverse.CalcGradR=true;
end


if CtrlVar.Inverse.ResetPersistentVariables
    ubP=[];
    vbP=[];
end

if ~isempty(ubP)
    F.ub=ubP;
    F.vb=vbP;
end


%% Reflection?


if contains(CtrlVar.Inverse.MinimisationMethod,"Ua")
    % pub and plb are enforced by the MATLAB optimization toolbox in a different way
    if CtrlVar.ReflectiveTransformation

        if ~isempty(pub)  && isempty(plb)

            iu=p>pub;
            p(iu)=pub(iu)+2*p(iu) ; % so if we had p(il)=plb(il) we get p(il)=plb(il)-2*p(il)=plb(il)

        elseif isempty(pub)  && ~isempty(plb)

            il=p>plb;
            p(il)=pub(il)+2*p(il) ; % so if we had p(il)=plb(il) we get p(il)=plb(il)-2*p(il)=plb(il)

        elseif ~isempty(pub)  && ~isempty(plb)

            %%
            % pub=[10 8]; plb=[1 2] ; p=[1 9] ;

            d=pub-plb;
            t=mod(p-plb,2*d);
            p=plb+min(t,2*d-t);
            %%

        end
    else
       
        p=kk_proj(p,pub,plb);
    end

end


%%




% The vector p contains the variables for which the inversion is being performed. So if the inversion is done over log(c)
% only, then p=log(C). And if the inversion is done over A, B and C then p=[A;B;C].

% Populate F with the current values in the vector p ahead of a call the the forward model.
F=p2F(CtrlVar,MUA,p,F,Meas,Priors);

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
[UserVar,RunInfo,F,l,dFduv]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);

if contains(CtrlVar.Inverse.Measurements,"-dhdt-")
    if isempty(F.as) || isempty(F.ab)
        [UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F);
    end
    [~,F.dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs) ;
end

%%
% The cost function, J), is split into a misfit (I) and a regularization term (R). These usually consist of further
% terms.
%
% Get the I and R terms, and the gradients if required.
if nargout==1
    R=Regularisation(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo) ;
    I=Misfit(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo,dFduv) ;
else
    [R,dRdp,ddRddp,RegOuts]=Regularisation(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo) ;
    [I,dIdp,ddIddp,MisfitOuts]=Misfit(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo,dFduv) ;
end




if nargout>1   % gradient needed
    dJdp=dRdp+dIdp;
end

if nargout>2  % Hessian needed
    if isempty(ddIddp)
        Hessian=ddRddp;
    else
        Hessian=ddRddp+ddIddp;
    end
end


if RunInfo.Forward.uvConverged
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


if nargout>3
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

end

