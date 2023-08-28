function [x,lambda,R2,r2,Slope0,dxNorm,dlambdaNorm,residual,g,h,output] = lsqLevenbergMarquardtUa(CtrlVar,fun,x,lambda,L,c)

%%
%
% Minimizes the norm of R where R is a vector subject the the
% constraints
%
%   L x = c
%
%  The function 'fun' should provide both R and the Jacobian, i.e.
%
%   [R,K]=fun(x) 
%
% where   K=grad J
%
% Note that R is a vector and K a matrix.  The value that is minimized is
%
%   R'*R
%
% Algorithm: Solves repeatedly the linearized min problem:
%
%  |R0+J dx|^2 
%
%
% Example:
%
%   lsqUaExample
%
%%
if isempty(CtrlVar) || ~isstruct(CtrlVar)

    ItMax=5;
    gTol=1e-20;
    dR2Tol=1e-3; 
    dxTol=1e-20;

    isLSQ=true;
    Normalize=false ;
    SaveIterate=true; 
    LevenbergMarquardt="auto" ; % "fixed"
    LMlambda=1 ;
    ScaleProblem=false;
    LMlambdaUpdateMethod=1; 


else

    ItMax=CtrlVar.lsqUa.ItMax ;
    gTol=CtrlVar.lsqUa.gTol ;
    dR2Tol=CtrlVar.lsqUa.dR2Tol ;
    dxTol=CtrlVar.lsqUa.dxTol ;

    isLSQ=CtrlVar.lsqUa.isLSQ ;
    Normalize=CtrlVar.lsqUa.Normalize;
    LevenbergMarquardt=CtrlVar.lsqUa.LevenbergMarquardt;
    LMlambda=CtrlVar.lsqUa.LMlambda0 ;
    ScaleProblem=CtrlVar.lsqUa.ScaleProblem; 

    LMlambdaUpdateMethod=CtrlVar.lsqUa.LMlambdaUpdateMethod ;

    SaveIterate=CtrlVar.lsqUa.SaveIterate;

end

R2Array=nan(ItMax+1,1) ;
g2Array=nan(ItMax+1,1) ;
dR2=[inf ; inf ] ; % stores the changes in R2=R'*R  over last two iterations

Slope0=nan ; % need to calculate final slope
nx=numel(x);

%% If contraints provided, make iterate feasable
if ~isempty(L) && ~isempty(c)   % if the user has not provided an initial estimate of lambda, but specifies constraints, set lambda=0

    BCres=norm(L*x-c);
    if BCres>1e-6   % make feasable
        x=L\c ;
    end
    if isempty(lambda)  || isnan(lambda)
        lambda=c*0;
    end
end


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
r2 = full(g'*g)/Normalisation;




g2Array(1)=r2;
R2Array(1)=R2;
if SaveIterate
    xVector=nan(numel(x),100);
    xVector(:,1)=x(:) ;
else
    xVector=[];
end
iteration=0 ;  

fprintf("\n\t Start lsqUa: \t  g=%g \t         r=%g \n \n",r2,R2)

while iteration <= ItMax

    iteration=iteration+1 ;

    if isLSQ
        KK=K'*K;
        H=2*KK;
        if LMlambda > 0
            if ScaleProblem
                D=sparse(1:nx,1:nx,spdiags(H,0));
            else
                D=speye(nx,nx);
            end
            H=H+D*LMlambda ;
        end
    else
        H=K;
    end

    [dx,dlambda]=solveKApe(H,L,g,h,x,lambda,CtrlVar);


   if isLSQ
       Q=2*R'*K*dx+dx'*KK*dx ;  % Quad approximation, based on unperturbed H
   else
       Q=R'*dx+dx'*H*dx/2 ;
   end
    
    x0=x ; lambda0=lambda ;

    x=x+dx ;
    lambda=lambda+dlambda ;



    R0=R ; K0=K ; g0=g ; h0=h ; % If I reject the step, I need to reuse these variables

    [R,K]=fun(x) ;

    if ~isempty(L)
        LTlambda=L'*lambda ;
        h =- (L*x-c);
    else
        LTlambda=0;
        h=[];
    end

    if isLSQ
        g =- (2*K'*R + LTlambda) ;
    else
        g =- (R + LTlambda) ;
    end


    g2Old=r2;  R2Old=R2;
    r2=full(g'*g)/Normalisation ;
    R2=full(R'*R);

    rho=(R2-R2Old)/Q;      % Actual reduction / Modelled Reduction
    g2Ratio=r2/g2Old ;
    dR2=[abs(R2-R2Old); dR2(1)] ; 


    if r2 > g2Old  % reject step
        x=x0 ; lambda=lambda0 ; R=R0 ; K=K0 ; g=g0 ; h=h0 ; r2=g2Old ; R2=R2Old ;
        
        if LMlambda==0
            LMlambda=1;
        else
            LMlambda=2*LMlambda ;
        end
        % fprintf(" step rejected \n")
        % fprintf(" %i",iteration)
        continue
    end


    dxNorm=norm(dx);
    dlambdaNorm=norm(dlambda);
    BCsNorm=norm(h) ;

    if isLSQ  && LevenbergMarquardt == "auto"
        if LMlambda==0
            LMlambda=1;
        end

        switch LMlambdaUpdateMethod
            % now update LMlambda based on agreement with Quad model
            case 1
                
                factor=max(1/3,1-(2*rho-1)^3) ;  % 
                LMlambda=factor*LMlambda ;

            case 2
                
                if rho>0.9  && rho<1.1       % Quad model very good, decrease lambda
                    LMlambda=LMlambda/10 ;
                elseif rho>0.75  && rho<1.5
                    LMlambda=LMlambda/2 ;
                elseif rho>0.1  || rho<10
                    LMlambda=1.5*LMlambda ;
                elseif rho < 0.1 || rho>10    % Quad model not good, increase lambda
                    LMlambda=2*LMlambda ;
                end

        end

    else
        LMlambda=nan;
    end

    g2Array(iteration+1)=r2;
    R2Array(iteration+1)=R2;
    if SaveIterate
        xVector(:,iteration+1)=x(:) ;
    end


    fprintf("lsqUa: \t it=%2i  \t     g0=%-13g \t     g1=%-13g \t         g1/g0=%-13g \t |R|^2=%-13g \t |dx|=%-13g \t |dl|=%-13g \t |BCs|=%-13g \t dr/Q=%-5f \t LMlambda=%g \n",iteration,g2Old,r2,g2Ratio,R2,dxNorm,dlambdaNorm,BCsNorm,rho,LMlambda)


    if r2 < gTol
        fprintf("lsqUa: Exiting iteration because |g|^2 within set tolerance of %g \n",gTol)
        break
    end

    if dxNorm < dxTol
        fprintf("lsqUa: Exiting iteration because change in x within the set tolerance of %g \n",dxTol)
        break
    end

    
    maxdR2=max(dR2);
    if maxdR2 < dR2Tol
        fprintf("lsqUa: Exiting iteration because max change in |R|^2 over last two iterations (%g) less than the set tolerance of %g \n",maxdR2,dR2Tol)
        break
    end
    
    if iteration >= ItMax
       fprintf("lsqUa: Exiting iteration because number of iterations has reached the set maximum of %i \n",ItMax)
       break

    end


end

fprintf("\n\t Exit lsqUa: \t  g=%g \t         r=%g \n \n",r2,R2)


residual=R ;
output.g2Array=g2Array;
output.R2Array=R2Array;
output.xVector=xVector; 
output.nIt=iteration; 



end

