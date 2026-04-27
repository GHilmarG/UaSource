



function [x,l,R2,r2,Slope0,dxNorm,dlNorm,residual,g,h,output] = lsqNewton(CtrlVar,fun,x,l,Aeq,beq)



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
isLSQ=true;
Normalize=false ;
SaveIterate=false;
lsqDogLeg="-Newton-Cauchy-";
CostMeasure="R2" ; % "r2"
StepString="  ";
InfoLevelNonLinIt= 1;

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

    if isfield(CtrlVar.lsqUa,"isLSQ")
        isLSQ=CtrlVar.lsqUa.isLSQ ;
    end

    if isfield(CtrlVar.lsqUa,"Normalize")
        Normalize=CtrlVar.lsqUa.Normalize;
    end
    if isfield(CtrlVar.lsqUa,"SaveIterate")
        SaveIterate=CtrlVar.lsqUa.SaveIterate;
    end

    if isfield(CtrlVar.lsqUa,"InfoLevelNonLinIt")
        InfoLevelNonLinIt=CtrlVar.InfoLevelNonLinIt;
    end

    if isfield(CtrlVar.lsqUa,"DogLeg")
        lsqDogLeg=CtrlVar.lsqUa.DogLeg ;
    end

    if isfield(CtrlVar.lsqUa,"CostMeasure")
        CostMeasure=CtrlVar.lsqUa.CostMeasure;
    end

    if isfield(CtrlVar,"CtrlVar.InfoLevelNonLinIt")
        InfoLevelNonLinIt=CtrlVar.InfoLevelNonLinIt;
    end

end


R2Array=nan(ItMax+1,1) ;  % value of the sum of squares of the FE right-hand system
r2Array=nan(ItMax+1,1) ;  % right-hand side of the lsq KKT system, first-order optimality
QArray=nan(ItMax+1,1) ;   % lsq quadratic approximation 
dxArray=nan(ItMax+1,1) ;
Slope0Array=nan(ItMax+1,1) ;
% WorkArray=nan(ItMax+1,1) ;
dR2=[inf ; inf ] ; % stores the changes in R2=R'*R  over last two iterations



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



g =- (2*K'*R + LTl) ;

if ~isempty(Aeq)
    h =- (Aeq*x-beq);
else
    h=[];
end

R2=full(R'*R);
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

while iteration <= ItMax

    iteration=iteration+1 ;

    K0=K ; R0=R; x0=x ; l0=l ;
    R20=R2;  r20=r2 ;  J0=J ;
    g0=g ; h0=h ; 
 
    KK0=2*(K0'*K0);

    gamma=nan; 
   

    % Solve the KKT Newton system. [dx;dl] is the Newton direction
    [dx,dl]=solveKApeSymmetric(KK0,Aeq,g0,h0,x0,l0,CtrlVar);

   [Qslope0,QgammaMin]=QgHlsq(R,K,Aeq,beq,x0,l0,dx,dl,gamma) ;

    CtrlVar.BacktrackIteration=iteration  ;

    
    CtrlVar.BacktrackingGammaMin=0.001;

    funcBackTrack=@(gamma) Jlsqfunc(CtrlVar,gamma,dx,dl,fun,Aeq,beq,x0,l0) ;

    
    J0=funcBackTrack(0) ; 
    J1=funcBackTrack(1) ; 

    CtrlVar.InfoLevelBackTrack=10000;  CtrlVar.InfoLevelNonLinIt=10 ;  CtrlVar.doplots=1 ;

 

    gammaNewton=1; 
    [gammamin,Jmin,BackTrackInfo]=BackTracking(Qslope0,gammaNewton,J0,J1,funcBackTrack,CtrlVar);

    x=x0+gammamin*dx ; l=l0+gammamin*dl ;
    
    [R,K]=fun(x) ;
    [Qslope0,~,Qreduction]=QgHlsq(R,K,Aeq,beq,x,l,dx,dl,gammamin) ; % this is the slope at the end of the line search in the direction of [dx;dl]

    
    if ~isempty(Aeq)
        h =- (Aeq*x-beq);
        g =- (2*K'*R + Aeq'*l) ;
    else
       g =- 2*K'*R  ;
        h=[];
    end

    R2=full(R'*R);
    r2 = full([g;h]'*[g;h])/Normalisation;

    
    r2Ratio=r2/r20 ;
    R2Ratio=R2/R20 ;
    dR2=[abs(R2-R20); dR2(1)] ;

    rho=(R20-R2)/Qreduction;

    dxNorm=norm(dx);
    dlNorm=norm(dl);
    BCsNorm=norm(h) ;


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

        fprintf("lsqUa: \t it=%2i%s  \t     |R|^2=%-13g \t     |R|^2/|R0|^2=%-13g \t gamma=%-13g \t |r|^2=%-13g \t |dx|=%-13g \t |dl|=%-13g \t |BCs|=%-13g \t dr/Q=%-5f \t slope0 =%g \n",...
            iteration,StepString,R2,R2Ratio,gammamin,r2,dxNorm,dlNorm,BCsNorm,rho,Qslope0)

    elseif CostMeasure=="r2"

        fprintf("lsqUa:%2i%s  \t     |r|^2=%-13g \t   |r0|^2=%-13g \t   |r|^2/|r0|^2=%-13g \t gamma=%-13g \t |R|^2=%-13g \t |dx|=%-13g \t |dl|=%-13g \t |BCs|=%-13g \t dr/Q=%-5f \t slope0 =%g \n",...
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

Slope0=full(2*R'*K*dx) ;
Slope0Array(iteration+1)=Slope0;  % This is the slope in the direction dx based on final R and K values

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
ylabel("$\|R\|^2$",Interpreter="latex")
lg=legend(Interpreter="latex");

%%



end



