



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
dxArray=nan(ItMax+1,1) ;
Slope0Array=nan(ItMax+1,1) ;
WorkArray=nan(ItMax+1,1) ;
dR2=[inf ; inf ] ; % stores the changes in R2=R'*R  over last two iterations

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


while iteration <= ItMax

    iteration=iteration+1 ;

    K0=K ; R0=R; x0=x ; lambda0=lambda ; h0=h ; g0=g;
    R20=R2;  g20=g2 ; 
    g2=nan ;

    if isLSQ
        KK0=K0'*K0;
        H0=2*KK0;

    else
        H0=K0;
        KK0=K0'*K0;
    end


    CtrlVar.BacktrackIteration=iteration  ;

    %% Newton Step, with possible backtracking
    CtrlVar.BacktracFigName="Newton";
    CtrlVar.BacktrackStepRatio=1e-2;              
    StepString="N ";
    [R2,x,lambda,dx,dlambda,Slope0,gammamin,BackTrackInfo,gammaEst,exitflag]=lsqStepUa(CtrlVar,fun,x0,lambda0,L,c,H0,R20,K0,R0,g0,h0,KK0) ;
    
    %%
    if exitflag ==1
        fprintf("lsqUa: Exiting iteration because slope at origin in Newton line search positive (Slope=%g) \n",Slope0)
        gammamin=0; 
        

    end

    %% Cauchy Step, with possible backtracking
    % Should I try Cauchy step?

    R2ratio=R2/R20;
    if R2ratio >0.9 || gammamin < 0.001
        TryCauchyStep=true;
    else
        TryCauchyStep=false;
    end

    if TryCauchyStep

        xN=x ; % dont' loose the Newton solution, in case a Cauchy-to-Newton path is required
        lambdaN=lambda ;
        R2N=R2; 

        I0=speye(nx) ;
        CtrlVar.BacktracFigName="Cauchy";
        CtrlVar.BacktrackStepRatio=1e-5;
        [R2C,xC,lambdaC,dxC,dlambdaC,Slope0C,gammaminC,BackTrackInfoC,gammaEst,exitflag]=lsqStepUa(CtrlVar,fun,x0,lambda0,L,c,I0,R20,K0,R0,g0,h0,KK0) ;
        if exitflag ==1

            fprintf("lsqUa: Exiting iteration because slope at origin in Cauchy line search positive (Slope=%g) \n",Slope0)
            break
            

        end

        if R2C < R2N
            fprintf("Cauchy step outperformes Newton. R2C/R2N=%g \n",R2C/R2N)
            StepString="C ";
            x=xC ; lambda=lambdaC  ; Slope0=Slope0C ; gammamin=gammaminC ;
            R2=R2C ;
            % Cauchy-to-Newton
            %% Cauchy-to-Newton search?
            % Should I try Cauchy step?

            R2ratio=R2/R20;
            if R2ratio >0.9 && gammaminC/gammaEst > 0.9  % Only to Caucy-to-Newton if not sufficient reduction, AND min found close to Cauchy Point
                TryCauchy2Newton=true;
            else
                TryCauchy2Newton=false;
            end

            if TryCauchy2Newton
                x0=xC ; lambda0=lambdaC ;  dx=xN-xC  ;  dlambda=lambdaN-lambdaC;

                funcCN=@(gamma) R2func(gamma,dx,dlambda,fun,x0,lambda0) ;
                TolX=0.01;
                options = optimset('Display','iter','TolX',TolX,'OutputFcn',@outfunFminbnd);
                [gammaminCN,R2CN,exitflag,outputfminbnd]=fminbnd(funcCN,0,1,options) ;

                if R2CN < R2

                    fprintf("Cauchy-to-Newton step outperformes. R2C/R2N=%g \n",R2C/R2N)

                    x=x0+gammaminCN*dx ;lambda=lambda0+gammaminCN*dlambda ;  R2=R2CN ; 

                    if CtrlVar.InfoLevelBackTrack>=1000

                        [stop,Outs]=outfunFminbnd();
                        PlotCauchy2NewtonPath(CtrlVar,Outs.x,Outs.f,R2C,R2N,gammaminCN,R2CN,R20);

                    end
                end
            end
        end
    end


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



