function [x,lambda,resnorm,residual,g,h,output] = lsqUa(CtrlVar,fun,x,lambda,L,c)


%
%   [H L' ]  [dx]  = [g]
%   [L  0 ]  [l]     [h]
%
%

if isempty(CtrlVar) || isnan(CtrlVar)  

    ItMax=5;
    tol=1e-6;
    isLSQ=true;
    LevenbergMarquardt=nan; 
    Normalize=false ;

else

    ItMax=CtrlVar.lsqUa.ItMax ;
    tol=CtrlVar.lsqUa.tol ;
    isLSQ=CtrlVar.lsqUa.isLSQ ;
    LevenbergMarquardt=CtrlVar.lsqUa.LevenbergMarquardt; 
    Normalize=CtrlVar.lsqUa.Normalize;

end

if ~isempty(L) && ~isempty(c)   % if the user has not provided an initial estimate of lambda, but specifies constraints, set lambda=0
    if isempty(lambda)  || isnan(lambda)
        lambda=c*0;
    end
end

if ~isempty(L)
    LTlambda=L'*lambda ;
else
    LTlambda=0;
end

% Normalisation
if Normalize
    x0=x*0 ;
    R=fun(x0) ;

    g =- (R + LTlambda) ;
    Normalisation=full(g'*g);
else
    Normalisation=1;

end


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



r = full(g'*g)/Normalisation;
r0=r ;



rVector=nan(100,1) ;
iteration=0 ;

while iteration < ItMax  && r > tol

    iteration=iteration+1 ;

    if isLSQ
        H=2*(K'*K) ;
        if  ~isnan(LevenbergMarquardt)
            H=H+speye(numel(x))*CtrlVar.lsqUa.LevenbergMarquardt ;
        end
    else
        H=K;
    end

    [dx,dlambda]=solveKApe(H,L,g,h,x,lambda,CtrlVar);

    x=x+dx ;
    lambda=lambda+dlambda ;

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
    rVector(iteration+1)=r;

    R2=full(R'*R);
    dxNorm=norm(dx)/norm(x);
    dlambdaNorm=norm(dlambda)/norm(lambda);
    BCsNorm=norm(h) ;

    fprintf("\t      lsqUa: \t it=%i  \t     r0=%-10g \t     r1=%-10g \t         r1/r0=%-10g \t R2=%-10g \t |dx|=%-10g \t |dl|=%-10g \t |BCs|=%g \n",iteration,r0,r,r/r0,R2,dxNorm,dlambdaNorm,BCsNorm)

end

fprintf("\n\t Exit lsqUa: \t it=%i  \t     r0=%g \t     r1=%g \t         r1/r0=%g \n \n",iteration,r0,r,r/r0)

resnorm=r;
residual=R ;
exitflag=0;
output.rVector=rVector;





end

