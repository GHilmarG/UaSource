



function [x,lambda,R2,r2,Slope0,dxNorm,dlambdaNorm,residual,g,h,output] = lsqDogLegUa(CtrlVar,fun,x,lambda,L,c)






ItMax=5;
gTol=1e-20;
dR2Tol=1e-3;
dxTol=1e-20;
isLSQ=true;
Normalize=false ;
SaveIterate=false;
lsqDogLeg="-Newton-Cauchy-";
CostMeasure="R2" ; % "r2"

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

end


R2Array=nan(ItMax+1,1) ;
r2Array=nan(ItMax+1,1) ;
JArray=nan(ItMax+1,1) ;
dxArray=nan(ItMax+1,1) ;
Slope0Array=nan(ItMax+1,1) ;
% WorkArray=nan(ItMax+1,1) ;
dR2=[inf ; inf ] ; % stores the changes in R2=R'*R  over last two iterations



%% If contraints provided, make iterate feasable

if ~isempty(L) && ~isempty(c)   % if the user has not provided an initial estimate of lambda, but specifies constraints, set lambda=0

    BCres=norm(L*x-c);
    if BCres>1e-6   % make feasable
        x=L\c ;
    end
    if isempty(lambda)  || anynan(lambda)
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

% fprintf("\n\t Start lsqUa: \t  r2=%g \t         R2=%g \n \n",r2,R2)


%%
if contains(lsqDogLeg,"-Cauchy-")
    TryCauchyStep=true;
else
    TryCauchyStep=false;
end

if contains(lsqDogLeg,"-Newton-")
    TryNewtonStep=true;
else
    TryNewtonStep=false;
end

while iteration <= ItMax

    iteration=iteration+1 ;

    K0=K ; R0=R; x0=x ; lambda0=lambda ; h0=h ; g0=g;
    R20=R2;  r20=r2 ;  J0=J ;
    r2=nan ; JminN=nan ; JminC=nan ; R2minCN=nan ; Slope0=nan ;

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
        CtrlVar.lsqUa.Step="-Newton-";
        CtrlVar.BacktracFigName="Newton";
        CtrlVar.BacktrackStepRatio=1e-2;
        CtrlVar.LineSearchAllowedToUseExtrapolation=true;



        [JminN,dxN,dlambdaN,gammaminN,Slope0N,BackTrackInfo,gammaEstN,exitflag]=lsqStepUa(CtrlVar,fun,x0,lambda0,K0,R0,L,c);

        % xN=x0+dxN ; lambdaN=lambda0+dlambdaN ;
        xN=x0+gammaminN*dxN ; lambdaN=lambda0+gammaminN*dlambdaN ;
        if JminN < J
            J=JminN ;
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
            if contains(lsqDogLeg,"-Cauchy-")
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



    Jratio=J/J0;

    if ~contains(lsqDogLeg,"-Cauchy-")
        TryCauchyStep=false;
    elseif ~TryCauchyStep  % Even though Newton step did result in reduction, I might still want to try the Cauchy step as well if the step size or the reduction was small
        if Jratio >0.9 || gammamin < 0.001
            TryCauchyStep=true;
        end
    end

    if TryCauchyStep

        CtrlVar.lsqUa.Step="-Cauchy-";
        
        CtrlVar.BacktracFigName="Cauchy";
        CtrlVar.BacktrackStepRatio=1e-5;
        CtrlVar.LineSearchAllowedToUseExtrapolation=true;
      
        
        [JminC,dxC,dlambdaC,gammaminC,Slope0C,BackTrackInfo,gammaEstC,exitflag]=lsqStepUa(CtrlVar,fun,x0,lambda0,K0,R0,L,c);
      
        xC=x0+gammaminC*dxC ; lambdaC=lambda0+gammaminC*dlambdaC ;
        if JminC < J
            J=JminC ;
            x=x0+gammaminC*dxC ;
            lambda=lambda0+gammaminC*dlambdaC ;
            Slope0=Slope0C ;  % Even if I then do the C2N step, this will be the estimate for Slop0, ie in the Cauchy direction.
            StepString="C ";
            TryCauchy2Newton=false;
            gammamin=gammaminC ;
            if contains(lsqDogLeg,"-Newton-")
                fprintf("Cauchy step outperformes Newton. R2minC/R2minN=%g \n",JminC/JminN  )
            end
        end
           

        if exitflag ==1

            fprintf("lsqUa: Exiting iteration because slope at origin in Cauchy line search positive (Slope=%g) \n",Slope0)

        end

        % Try Cauchy-to-Newton?
        Jratio=J/J0;
        % Only to C2N if not sufficient reduction alread, both Newton and Cauchy have been performed and the min in the Cauchy step
        % was found close to the min of the Quad approximation
        if Jratio>0.9 &&  ~isnan(JminN) && ~isnan(JminC) && gammaminC/gammaEstC > 0.9
            TryCauchy2Newton=true;
        else
            TryCauchy2Newton=false;
        end

        if TryCauchy2Newton

            dxC2N=xN-xC  ;  dlambdaC2N=lambdaN-lambdaC;

            % funcCN=@(gamma) R2func(gamma,dxC2N,dlambdaC2N,fun,xC,lambdaC) ;


            funcCN=@(gamma)  Jlsqfunc(CtrlVar,gamma,dxC2N,dlambdaC2N,fun,L,c,xC,lambdaC);


            TolX=0.01;
            if InfoLevelNonLinIt >= 10
                options = optimset('Display','iter','TolX',TolX,'OutputFcn',@outfunFminbnd);
            else
                options = optimset('Display','off','TolX',TolX,'OutputFcn',@outfunFminbnd);
            end


            [gammaminCN,R2CN,exitflag,outputfminbnd]=fminbnd(funcCN,0,1,options) ;
           
            if CtrlVar.InfoLevelBackTrack>=1000

                [stop,Outs]=outfunFminbnd();
                PlotCauchy2NewtonPath(CtrlVar,Outs.x,Outs.f,JminC,JminN,gammaminCN,R2CN,R20);

            end

            if R2CN < R2

                fprintf("Cauchy-to-Newton step outperformes. R2minCN/R2=%g \n",JminC/R2)

                gammamin=gammaminCN ;
                x=xC+gammaminCN*dxC2N ;lambda=lambdaC+gammaminCN*dlambdaC2N ;


            end
        end
    end
    
    dx=x-x0 ; dlambda=lambda-lambda0 ;



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

    r2=full([g;h]'*[g;h])/Normalisation ;


    Q=2*R0'*K0*dx+dx'*KK0*dx ;  % Quad approximation, based on unperturbed H
    if isLSQ
        rho=(R2-R20)/Q;      % Actual reduction / Modelled Reduction
    else
        % if isempty(L)
        %     Q=R0'*dx + dx'*H0*dx/2 ;
        % else
        %     Q=(R0+L'*lambda0)'*dx/2+(L*x0-c)'*dlambda/2 ...
        %         + dx'*(H0*dx+L'*dlambda)/2 + dlambda'*L*dx/2;
        % end
        rho=(r2-r20)/Q;      % Actual reduction / Modelled Reduction
    end



%%

    
    r2Ratio=r2/r20 ;
    R2Ratio=R2/R20 ;
    dR2=[abs(R2-R20); dR2(1)] ;


    dxNorm=norm(dx);
    dlambdaNorm=norm(dlambda);
    BCsNorm=norm(h) ;


    r2Array(iteration+1)=r2;
    R2Array(iteration+1)=R2;
    dxArray(iteration)=dxNorm ;
    Slope0Array(iteration)=Slope0;   % This is the slope based on R0, K0 and dx. Note the slope in dx direction at the end of the step
    % If doing a line search, the slope at the end of the step should always be close to zero in
    % the direction dx.


  %  WorkArray(iteration+1)=[dx;dlambda]'*[g ; h] ;

    if SaveIterate
        xVector(:,iteration+1)=x(:) ;
    end


    if CostMeasure=="R2"

        fprintf("lsqUa: \t it=%2i%s  \t     |R|^2=%-13g \t     |R|^2/|R0|^2=%-13g \t gamma=%-13g \t |r|^2=%-13g \t |dx|=%-13g \t |dl|=%-13g \t |BCs|=%-13g \t dr/Q=%-5f \t slope0 =%g \n",...
            iteration,StepString,R2,R2Ratio,gammamin,r2,dxNorm,dlambdaNorm,BCsNorm,rho,Slope0)

    elseif CostMeasure=="r2"

        fprintf("lsqUa:%2i%s  \t     |r|^2=%-13g \t   |r0|^2=%-13g \t   |r|^2/|r0|^2=%-13g \t gamma=%-13g \t |R|^2=%-13g \t |dx|=%-13g \t |dl|=%-13g \t |BCs|=%-13g \t dr/Q=%-5f \t slope0 =%g \n",...
            iteration,StepString,r2,r20,r2Ratio,gammamin,R2,dxNorm,dlambdaNorm,BCsNorm,rho,Slope0)

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
    residual=R ;
else
    residual=g;
end

output.r2Array=r2Array;
output.R2Array=R2Array;
output.dxArray=dxArray;
output.Slope0Array=Slope0Array;
% output.WorkArray=WorkArray;

output.xVector=xVector;
output.nIt=iteration;

end



