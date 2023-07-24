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
    LevenbergMarquardt=nan; 
    Normalize=false ;
    SaveIterate=true; 

else

    ItMax=CtrlVar.lsqUa.ItMax ;
    tol=CtrlVar.lsqUa.tol ;
    isLSQ=CtrlVar.lsqUa.isLSQ ;
    LevenbergMarquardt=CtrlVar.lsqUa.LevenbergMarquardt; 
    Normalize=CtrlVar.lsqUa.Normalize;
    SaveIterate=CtrlVar.lsqUa.SaveIterate;

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
rVector(1)=r;
if SaveIterate
    xVector=nan(numel(x),100);
    xVector(:,1)=x(:) ;
else
    xVector=[];
end
iteration=0 ;

while iteration < ItMax  && r > tol

    iteration=iteration+1 ;

    if isLSQ
        H=2*(K'*K) ;
        if  ~isnan(LevenbergMarquardt)
            H=H+speye(numel(x))*LevenbergMarquardt ;
        end
    else
        H=K;
    end

    [dx,dlambda]=solveKApe(H,L,g,h,x,lambda,CtrlVar);

    x=x+dx ;
    lambda=lambda+dlambda ;

    if SaveIterate
        xVector(:,iteration+1)=x(:) ;
    end

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

    fprintf("\t      lsqUa: \t it=%2i  \t     r0=%-13g \t     r1=%-13g \t         r1/r0=%-13g \t R2=%-13g \t |dx|=%-13g \t |dl|=%-13g \t |BCs|=%g \n",iteration,r0,r,r/r0,R2,dxNorm,dlambdaNorm,BCsNorm)

end

fprintf("\n\t Exit lsqUa: \t it=%i  \t     r0=%g \t     r1=%g \t         r1/r0=%g \n \n",iteration,r0,r,r/r0)

resnorm=r;
residual=R ;
exitflag=0;
output.rVector=rVector;
output.xVector=xVector; 
output.nIt=iteration; 



end

