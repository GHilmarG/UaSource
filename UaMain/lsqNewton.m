



function [x,l,R2,r2,Qslope0,dxNorm,dlNorm,residual,g,h,output] = lsqNewton(CtrlVar,fun,x,l,Aeq,beq)



%%
%
% The function fun returns
%
% [R,K,outs]=fun(x)
%
% where R is a vector valued function, and the objective function if R'*R  
%
% K is the Jacobian of R
%
%%


ItMax=5;
gTol=1e-20;
dR2Tol=1e-3;
dxTol=1e-20;
Normalize=false ;
SaveIterate=false;
CostMeasure="R2" ; % "r2"
StepString="  ";

if ~isempty(CtrlVar) && isstruct(CtrlVar) && isfield(CtrlVar,"lsqUa")

    if isfield(CtrlVar.lsqUa,"ItMax")
        ItMax=CtrlVar.lsqUa.ItMax ;
    end

    if isfield(CtrlVar.lsqUa,"gTol")
        gTol=CtrlVar.lsqUa.gTol ;
    end

    if isfield(CtrlVar.lsqUa,"dR2Tol")

        dR2Tol=CtrlVar.lsqUa.dR2Tol ;
    end

    if isfield(CtrlVar.lsqUa,"dxTol")
        dxTol=CtrlVar.lsqUa.dxTol ;
    end

  
    if isfield(CtrlVar.lsqUa,"Normalize")
        Normalize=CtrlVar.lsqUa.Normalize;
    end
    if isfield(CtrlVar.lsqUa,"SaveIterate")
        SaveIterate=CtrlVar.lsqUa.SaveIterate;
    end

 

    if isfield(CtrlVar.lsqUa,"CostMeasure")
        CostMeasure=CtrlVar.lsqUa.CostMeasure;
    end

end


R2Array=nan(ItMax+1,1) ;  % value of the sum of squares of the FE right-hand system
r2Array=nan(ItMax+1,1) ;  % right-hand side of the lsq KKT system, first-order optimality
QArray=nan(ItMax+1,1) ;   % lsq quadratic approximation 
dxArray=nan(ItMax+1,1) ;
Slope0Array=nan(ItMax+1,1) ;
% WorkArray=nan(ItMax+1,1) ;
dR2=[inf ; inf ] ; % stores the changes in R2=0.5*R'*R  over last two iterations



%% If constraints provided, make iterate feasible

if ~isempty(Aeq) && ~isempty(beq)   % if the user has not provided an initial estimate of l, but specifies constraints, set l=0

    BCres=norm(Aeq*x-beq);
    if BCres>1e-6   % make feasible
        x=Aeq\beq ;
    end
    if isempty(l)  || anynan(l)
        l=beq*0;
    end
end

if ~isempty(Aeq)
    LTl=Aeq'*l ;
else
    LTl=0;
end
%%

%% Normalisation
if Normalize
    x0=x*0 ;
    R=fun(x0) ;

    g =- (R + LTl) ;
    Normalisation=full(g'*g);
else
    Normalisation=1;

end

% Evaluate cost function, don't solve system
[R,K]=fun(x) ;



g =- (K'*R + LTl) ;

if ~isempty(Aeq)
    h =- (Aeq*x-beq);
else
    h=[];
end

R2=0.5*full(R'*R);
r2 = full([g;h]'*[g;h])/Normalisation;

if CostMeasure=="R2"
    J=R2;
elseif CostMeasure=="r2"
    J=r2;
else
    error("what cost measure?")
end

r2Array(1)=r2;
R2Array(1)=R2;


if SaveIterate
    xVector=nan(numel(x),100);
    xVector(:,1)=x(:) ;
else
    xVector=[];
end
iteration=0 ;

%

Delta=nan ; % trust-region radius (if used)
alpha=0;
TrustRegionSubProblem=true;



useMATLABquadprog=true;

while iteration <= ItMax

    iteration=iteration+1 ;

    x0=x ;
    l0=l ;

    R20=R2;
    r20=r2 ;


    KK=K'*K;

    gamma=nan;



    if TrustRegionSubProblem
        if ~isnan(Delta)
            E=speye(size(KK))  ;
            alpha=TrustRegionSubproblem(KK,E,g,alpha,Delta);
            E=speye(size(KK)) ;
            KK=KK+alpha*E;
        end
    end




    if useMATLABquadprog

        options=optimoptions('quadprog',display='iter') ;

        [dx,fval,exitflag,output,lambda]=quadprog(KK,-g,[],[],Aeq,h,[],[],[],options);
        dl=lambda.eqlin;

    else

        % Solve the KKT Newton system. [dx;dl] is the Newton direction

        [dx,dl]=solveKApeSymmetric(KK,Aeq,g,h,x0,l0,CtrlVar);

        dl=full(dl);
    end


  
    if isnan(Delta) || alpha==0 
        Delta=norm(dx); % If I have not set the trust-region radius, set it to the norm of the previous solution, which here will be the full Newton step
    end

    [Qslope0,QgammaMin,Qreduction]=QgHlsq(R,K,Aeq,beq,x0,l0,dx,dl,gamma) ;


    if Qslope0>0
    
        fprintf("Slope of quadradic model at origin in the search direction is positive, with slope0=%g \n",Qslope0)

    end

    CtrlVar.BacktrackIteration=iteration  ;


    CtrlVar.BacktrackingGammaMin=0.001;

    funcBackTrack=@(gamma) Jlsqfunc(CtrlVar,gamma,dx,dl,fun,Aeq,beq,x0,l0) ;


    J0=funcBackTrack(0) ; 
    J1=funcBackTrack(1) ; 

    CtrlVar.InfoLevelBackTrack=100;  CtrlVar.InfoLevelNonLinIt=10 ;  CtrlVar.doplots=1 ;



    gammaNewton=1;
    CtrlVar.LineSearchAllowedToUseExtrapolation=true;
    [gammamin,Jmin,BackTrackInfo]=BackTracking(Qslope0,gammaNewton,J0,J1,funcBackTrack,CtrlVar);

    if isnan(gammamin)
        gammamin=0;  % this can happen if, for example, the slope at gamma=0 is positive.
        % this might also be a numerical issue if the slope is very small
        % Basically, either exit the iteration, or try new Delta/alpha update when using the Trust-Region algorithm
        rho=0;
        if TrustRegionSubProblem
            Delta=TrustRegionRadiusUpdate(Delta,rho) ;
            fprintf("Backtracking resulted in gamma=%g . Decreasing trust-region radius to Delta=%g and repeating step. \n",gammamin,Delta)

            r2Array(iteration+1)=r2;
            R2Array(iteration+1)=R2;
            dxArray(iteration)=dxNorm ;

            continue
        else
            fprintf("Backtracking resulted in gamma=%g . Exiting Newton iteration. \n",gammamin)
            break
        end

    end


    [~,~,Qreduction]=QgHlsq(R,K,Aeq,beq,x0,l0,dx,dl,gammamin) ;

    x=x0+gammamin*dx ; l=l0+gammamin*dl ;
    
    % Now I've updated x and I need to recalculate K and R
    [R,K]=fun(x) ; % This is the only call within the NR loop.
  
    


    if ~isempty(Aeq)
        h =- (Aeq*x-beq);
        g =- (K'*R + Aeq'*l) ;
    else
        g =- K'*R  ;
        h=[];
    end

    dxNorm=norm(dx);
    dlNorm=norm(dl);
    BCsNorm=norm(h) ;

    R2=0.5*full(R'*R);
    r2 = 0.5*full([g;h]'*[g;h])/Normalisation;


    r2Ratio=r2/r20 ;
    R2Ratio=R2/R20 ;
    dR2=[abs(R2-R20); dR2(1)] ;


    if gammamin==0 || Qreduction==0
        rho=0;
    else
        rho=(R2-R20)/Qreduction;
    end

  



    r2Array(iteration+1)=r2;
    R2Array(iteration+1)=R2;
    dxArray(iteration)=dxNorm ;
    Slope0Array(iteration)=Qslope0;   % This is the slope based on R0, K0 and dx. Note the slope in dx direction at the end of the step
    % If doing a line search, the slope at the end of the step should always be close to zero in
    % the direction dx.


  %  WorkArray(iteration+1)=[dx;dl]'*[g ; h] ;

    if SaveIterate
        xVector(:,iteration+1)=x(:) ;
    end


    if CostMeasure=="R2"

        fprintf("lsqUa: \t it=%2i%s  \t     |R|^2=%-13g \t     |R|^2/|R0|^2=%-13g \t gamma=%-13g \t |r|^2=%-13g \t |dx|=%-13g \t |dl|=%-13g \t |BCs|=%-13g \t rho=%-5f \t slope0 =%g \n",...
            iteration,StepString,R2,R2Ratio,gammamin,r2,dxNorm,dlNorm,BCsNorm,rho,Qslope0)

    elseif CostMeasure=="r2"

        fprintf("lsqUa:%2i%s  \t     |r|^2=%-13g \t   |r0|^2=%-13g \t   |r|^2/|r0|^2=%-13g \t gamma=%-13g \t |R|^2=%-13g \t |dx|=%-13g \t |dl|=%-13g \t |BCs|=%-13g \t rho=%-5f \t slope0 =%g \n",...
            iteration,StepString,r2,r20,r2Ratio,gammamin,R2,dxNorm,dlNorm,BCsNorm,rho,Qslope0)

    else

        error("what cost measure?")

    end

    if r2 < gTol
        fprintf("lsqUa: Exiting iteration because |g|^2=%g within set tolerance of %g \n",r2,gTol)
        break
    end

    if dxNorm < dxTol
        fprintf("lsqUa: Exiting iteration because change in |x|=%g within the set tolerance of %g \n",dxNorm,dxTol)
        break
    end


    maxdR2=max(dR2);
    if maxdR2 < dR2Tol
        fprintf("lsqUa: Exiting iteration because max change in |R|^2=%g over last two iterations, less than the set tolerance of %g \n",maxdR2,dR2Tol)
        break
    end

    if iteration >= ItMax
        fprintf("lsqUa: [\b Exiting]\b  iteration because number of iterations has reached the set maximum of %i \n",ItMax)
        break

    end


end


Slope0Array(iteration+1)=Qslope0;  % This is the slope in the direction dx based on final R and K values

% fprintf("\n\t Exit lsqUa: \t  |g|^2=%g \t    slope=%g \t     |R|^2=%g \n \n",r2,Slope0,R2)

if CostMeasure=="R2"
    residual=R2;
elseif CostMeasure=="r2"
    residual=r2;
else
    error("what cost measure?")
end


output.r2Array=r2Array;
output.R2Array=R2Array;
output.dxArray=dxArray;
output.Slope0Array=Slope0Array;
% output.WorkArray=WorkArray;

output.xVector=xVector;
output.nIt=iteration;
%output.fun=funOuts ; % missing this, need to add


%%
FigNL=FindOrCreateFigure("Non-lin Convergence") ;  clf(FigNL)
hold off
yyaxis left
plot(0:iteration,output.r2Array(1:iteration+1),"bo-",DisplayName="$r^2$ (first-order optimality)")  
FigNL.CurrentAxes.YScale="log"   ;
ylabel("$r^2$, first-order optimality",Interpreter="latex")
hold on 
yyaxis right
plot(0:iteration,output.R2Array(1:iteration+1),"ro-",DisplayName="$\|R\|^2$")  
FigNL.CurrentAxes.YScale="log"   ;
ylabel("$\|R\|^2$",Interpreter="latex")
lg=legend(Interpreter="latex");

%%



end



