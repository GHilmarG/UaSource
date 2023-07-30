



function [x,lambda,R2,Slope0,dxNorm,dlambdaNorm,g2,residual,g,h,output] = lsqDogLegUa(CtrlVar,fun,x,lambda,L,c)




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
    InfoLevelNonLinIt=1;

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
    InfoLevelNonLinIt=CtrlVar.InfoLevelNonLinIt;

end

R2Array=nan(ItMax+1,1) ;
g2Array=nan(ItMax+1,1) ;
dxArray=nan(ItMax+1,1) ;
Slope0Array=nan(ItMax+1,1) ;
WorkArray=nan(ItMax+1,1) ;
dR2=[inf ; inf ] ; % stores the changes in R2=R'*R  over last two iterations

nx=numel(x); nlambda=numel(lambda) ;

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
%%

%% Normalisation
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

R2=full(R'*R);
g2 = full(g'*g)/Normalisation;


g2Array(1)=g2;
R2Array(1)=R2;

if SaveIterate
    xVector=nan(numel(x),100);
    xVector(:,1)=x(:) ;
else
    xVector=[];
end
iteration=0 ;

fprintf("\n\t Start lsqUa: \t  g=%g \t         r=%g \n \n",g2,R2)


%%
if contains(CtrlVar.lsqDogLeg,"-Cauchy-")
    TryCauchyStep=true;
else
    TryCauchyStep=false;
end

if contains(CtrlVar.lsqDogLeg,"-Newton-")
    TryNewtonStep=true;
else
    TryNewtonStep=false;
end

while iteration <= ItMax

    iteration=iteration+1 ;

    K0=K ; R0=R; x0=x ; lambda0=lambda ; h0=h ; g0=g;
    R20=R2;  g20=g2 ; 
    g2=nan ; R2minN=nan ; R2minC=nan ; R2minCN=nan ; Slope0=nan ; 

    if isLSQ
        KK0=K0'*K0;
        H0=2*KK0;

    else
        H0=K0;
        KK0=K0'*K0;
    end


    CtrlVar.BacktrackIteration=iteration  ;

    if TryNewtonStep

        %% Newton Step, with possible backtracking
        CtrlVar.BacktracFigName="Newton";
        CtrlVar.BacktrackStepRatio=1e-2;
        CtrlVar.LineSearchAllowedToUseExtrapolation=false;

        [R2minN,dxN,dlambdaN,gammaminN,Slope0N,BackTrackInfo,gammaEstN,exitflag]=lsqStepUa(CtrlVar,fun,x0,lambda0,L,H0,R20,K0,R0,g0,h0,KK0) ;
        xN=x0+dxN ; lambdaN=lambda0+dlambdaN ;
        if R2minN < R2
            R2=R2minN ;
            x=x0+gammaminN*dxN ;
            lambda=lambda0+gammaminN*dlambdaN ;
            StepString="N ";
            TryCauchyStep=false;
            gammamin=gammaminN ;
            Slope0=Slope0N ;
        else
            x=x0; lambda=lambda0;
            StepString="  ";
            gammamin=0 ;
            if contains(CtrlVar.lsqDogLeg,"-Cauchy-")
                TryCauchyStep=true;
            else
                TryCauchyStep=false;
            end
        end

        if exitflag ==1  % do I need this?
            fprintf("lsqUa: Slope at origin in Newton line search positive (Slope=%g) \n",Slope0)
            gammamin=0;
        end

    end

    %% Cauchy Step, with possible backtracking
    % Should I try Cauchy step?



    R2ratio=R2/R20;
    if ~TryCauchyStep  % Even though Newton step did result in reduction, I might still want to try the Cauchy step as well if the step size or the reduction was small
        if R2ratio >0.9 || gammamin < 0.001
            TryCauchyStep=true;
        end
    end

    if TryCauchyStep

        I0=speye(nx) ;
        CtrlVar.BacktracFigName="Cauchy";
        CtrlVar.BacktrackStepRatio=1e-5;
        CtrlVar.LineSearchAllowedToUseExtrapolation=true;
        [R2minC,dxC,dlambdaC,gammaminC,Slope0C,BackTrackInfo,gammaEstC,exitflag]=lsqStepUa(CtrlVar,fun,x0,lambda0,L,I0,R20,K0,R0,g0,h0,KK0) ;
        xC=x0+dxC ; lambdaC=lambda0+dlambdaC ;
        if R2minC < R2
            R2=R2minC ;
            x=x0+gammaminC*dxC ;
            lambda=lambda0+gammaminC*dlambdaC ;
            Slope0=Slope0C ;  % Even if I then do the C2N step, this will be the estimate for Slop0, ie in the Cauchy direction.
            StepString="C ";
            TryCauchy2Newton=false;
            gammamin=gammaminC ;
            if contains(CtrlVar.lsqDogLeg,"-Newton-")
                fprintf("Cauchy step outperformes Newton. R2minC/R2minN=%g \n",R2minC/R2minN  )
            end
        else
            x=x0 ; lambda=lambda0 ;
            StepString="  ";
            gammamin=0 ;
            TryCauchy2Newton=true;
        end

        if exitflag ==1

            fprintf("lsqUa: Exiting iteration because slope at origin in Cauchy line search positive (Slope=%g) \n",Slope0)

        end

        % Try Cauchy-to-Newton?
        R2ratio=R2/R20;
        % Only to C2N if not sufficient reduction alread, both Newton and Cauchy have been performed and the min in the Cauchy step
        % was found close to the min of the Quad approximation
        if R2ratio>0.9 &&  ~isnan(R2minN) && ~isnan(R2minC) && gammaminC/gammaEstC > 0.9
            TryCauchy2Newton=true;
        end

        if TryCauchy2Newton

            dxC2N=xN-xC  ;  dlambdaC2N=lambdaN-lambdaC;

            funcCN=@(gamma) R2func(gamma,dxC2N,dlambdaC2N,fun,xC,lambdaC) ;
            TolX=0.01;
            if InfoLevelNonLinIt >= 10
                options = optimset('Display','iter','TolX',TolX,'OutputFcn',@outfunFminbnd);
            else
                options = optimset('Display','off','TolX',TolX,'OutputFcn',@outfunFminbnd);
            end


            [gammaminCN,R2CN,exitflag,outputfminbnd]=fminbnd(funcCN,0,1,options) ;

            if R2CN < R2

                fprintf("Cauchy-to-Newton step outperformes. R2minCN/R2=%g \n",R2minC/R2)

                gammamin=gammaminCN ;
                x=xC+gammaminCN*dxC2N ;lambda=lambdaC+gammaminCN*dlambdaC2N ;

                if CtrlVar.InfoLevelBackTrack>=1000

                    [stop,Outs]=outfunFminbnd();
                    PlotCauchy2NewtonPath(CtrlVar,Outs.x,Outs.f,R2minC,R2minN,gammaminCN,R2CN,R20);

                end
            end
        end
    end
    
    dx=x-x0 ; dlambda=lambda-lambda0 ; 

    if isLSQ
        Q=2*R0'*K0*dx+dx'*KK0*dx ;  % Quad approximation, based on unperturbed H
    else
        Q=R0'*dx+dx'*H0*dx/2 ;
    end


    [R,K]=fun(x) ;
    R2=full(R'*R);

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

    g2=full(g'*g)/Normalisation ;
    
    

   %%

    rho=(R2-R20)/Q;      % Actual reduction / Modelled Reduction
    g2Ratio=g2/g20 ;
    R2Ratio=R2/R20 ;
    dR2=[abs(R2-R20); dR2(1)] ;


    dxNorm=norm(dx);
    dlambdaNorm=norm(dlambda);
    BCsNorm=norm(h) ;


    g2Array(iteration+1)=g2;
    R2Array(iteration+1)=R2;
    dxArray(iteration)=dxNorm ; 
    Slope0Array(iteration)=Slope0;   % This is the slope based on R0, K0 and dx. Note the slope in dx direction at the end of the step 
                                     % If doing a line search, the slope at the end of the step should always be close to zero in
                                     % the direction dx.


    WorkArray(iteration+1)=[dx;dlambda]'*[g ; h] ;

    if SaveIterate
        xVector(:,iteration+1)=x(:) ;
    end


    fprintf("lsqUa: \t it=%2i%s  \t     |R|^2=%-13g \t     |R|^2/|R0|^2=%-13g \t gamma=%-13g \t |g|^2=%-13g \t |dx|=%-13g \t |dl|=%-13g \t |BCs|=%-13g \t dr/Q=%-5f \t slope0 =%g \n",...
        iteration,StepString,R2,R2Ratio,gammamin,g2,dxNorm,dlambdaNorm,BCsNorm,rho,Slope0)


    if g2 < gTol
        fprintf("lsqUa: Exiting iteration because |g|^2 within set tolerance of %g \n",gTol)
        break
    end

    if dxNorm < dxTol
        fprintf("lsqUa: Exiting iteration because change in |x|=%g within the set tolerance of %g \n",dxNorm,dxTol)
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

Slope0=2*R'*K*dx ;  
Slope0Array(iteration+1)=Slope0;  % This is the slope in the direction dx based on final R and K values

fprintf("\n\t Exit lsqUa: \t  |g|^2=%g \t    slope=%g \t     |R|^2=%g \n \n",g2,Slope0,R2)



residual=R ;
output.g2Array=g2Array;
output.R2Array=R2Array;
output.dxArray=dxArray;
output.Slope0Array=Slope0Array;
output.WorkArray=WorkArray;

output.xVector=xVector;
output.nIt=iteration;

end



