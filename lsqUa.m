function [x,lambda,resnorm,residual,g,h,output] = lsqUa(CtrlVar,fun,x,lambda,L,c)

%%
%
%
%
%   [H L' ]  [dx]  = [g]
%   [L  0 ]  [l]     [h]
%
%
% Example:
%
%   lsqUaExample
%
%%
if isempty(CtrlVar) || ~isstruct(CtrlVar)

    ItMax=5;
    tol=1e-6;
    isLSQ=true;
 
    Normalize=false ;
    SaveIterate=true; 
    LevenbergMarquardt="auto" ; % "fixed"
    LMlambda=1 ;

else

    ItMax=CtrlVar.lsqUa.ItMax ;
    tol=CtrlVar.lsqUa.tol ;
    isLSQ=CtrlVar.lsqUa.isLSQ ;
    
    Normalize=CtrlVar.lsqUa.Normalize;
    LevenbergMarquardt=CtrlVar.lsqUa.LevenbergMarquardt;
    LMlambda=CtrlVar.lsqUa.LMlambda0 ;

    SaveIterate=CtrlVar.lsqUa.SaveIterate;

end

rVector=nan(ItMax+1,1) ;
gVector=nan(ItMax+1,1) ;

if ~isempty(L) && ~isempty(c)   % if the user has not provided an initial estimate of lambda, but specifies constraints, set lambda=0
    if isempty(lambda)  || isnan(lambda)
        lambda=c*0;
    end
end


%% is feasable?

BCres=norm(L*x-c);

if BCres>1e-6   % make feasable
    x=L\c ;
end

%%

if ~isempty(L)
    LTlambda=L'*lambda ;
else
    LTlambda=0;
end

%% Normalisation
if Normalize
    x0=x*0 ;
    R=fun(x0) ;

    g =- (R + LTlambda) ;
    Normalisation=full(g'*g);
else
    Normalisation=1;

end

%%

% Evaluate cost function, don't solve system
[R,K]=fun(x) ;

if isLSQ
    g =- (2*K'*R + LTlambda) ;
else
    g =- (R + LTlambda) ;
end

if ~isempty(L)
    h =- (L*x-c);
else
    h=[];
end

R2=full(R'*R); 
r = full(g'*g)/Normalisation;
r0=r ;



gVector(1)=r;
rVector(1)=R2;
if SaveIterate
    xVector=nan(numel(x),100);
    xVector(:,1)=x(:) ;
else
    xVector=[];
end
iteration=0 ;  

fprintf("\n\t Start lsqUa: \t  g1=%g \t         r=%g \n \n",r,R2)

while iteration < ItMax  && r > tol

    iteration=iteration+1 ;

    if isLSQ

        H=2*(K'*K)+speye(numel(x))*LMlambda ;

    else
        H=K;
    end

    [dx,dlambda]=solveKApe(H,L,g,h,x,lambda,CtrlVar);

    % expected value based on linearisation of R
    
   % gEst= 2*K' * (R + K*dx) + L'*(lambda+dlambda) ; % should be zero
   % because this is the system I solve.
   % rEst=(gEst'*gEst)/Normalisation ;

   if isLSQ
       Q=2*R'*K*dx+dx'*(K'*K)*dx ;
   else
       Q=R'*dx+dx'*K*dx ;
   end
    
    x0=x ; lambda0=lambda ;

    x=x+dx ;
    lambda=lambda+dlambda ;



    R0=R ; K0=K ; g0=g ; h0=h ; 
    [R,K]=fun(x) ;

    if ~isempty(L)
        LTlambda=L'*lambda ;
    else
        LTlambda=0;
    end

    if isLSQ
        g =- (2*K'*R + LTlambda) ;
    else
        g =- (R + LTlambda) ;
    end

    if ~isempty(L)
        h =- (L*x-c);
    else
        h=[];
    end

    r0=r;
    r=full(g'*g)/Normalisation ;
    R20=R2;
    R2=full(R'*R);
    rho=(R2-R20)/Q;


    if isLSQ  && LevenbergMarquardt == "auto"

        if r > r0  % reject step
            x=x0 ; lambda=lambda0 ; R=R0 ; K=K0 ; g=g0 ; h=h0 ; dx=dx*0; r=r0 ; dlambda=dlambda*0 ;
            LMlambda=10*LMlambda ;
            if LMlambda==0
                LMlambda=1;
            end
            %fprintf(" step rejected \n")

        else

            if rho>0.9  && rho<1.1
                LMlambda=LMlambda/10 ;
            elseif rho>0.75  && rho<1.5
                LMlambda=LMlambda/2 ;
            elseif rho>0.1  || rho<10
                LMlambda=1.5*LMlambda ;
            elseif rho < 0.1 || rho>10
                LMlambda=2*LMlambda ;
            end
        end
    end

    gVector(iteration+1)=r;
    rVector(iteration+1)=R2;
    if SaveIterate
        xVector(:,iteration+1)=x(:) ;
    end

    R2=full(R'*R);
    dxNorm=norm(dx)/norm(x);
    dlambdaNorm=norm(dlambda)/norm(lambda);
    BCsNorm=norm(h) ;





    fprintf("\t      lsqUa: \t it=%2i  \t     g0=%-13g \t     g1=%-13g \t         g1/g0=%-13g \t r=%-13g \t |dx|=%-13g \t |dl|=%-13g \t |BCs|=%-13g \t dr/Q=%-5f \t LMlambda=%g \n",iteration,r0,r,r/r0,R2,dxNorm,dlambdaNorm,BCsNorm,rho,LMlambda)

end

fprintf("\n\t Exit lsqUa: \t  g1=%g \t         r=%g \n \n",r,R2)

resnorm=r;
residual=R ;
exitflag=0;
output.gVector=gVector;
output.rVector=rVector;
output.xVector=xVector; 
output.nIt=iteration; 



end

