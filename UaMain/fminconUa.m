function [x,J,exitflag,output] = fminconUa(func,x0,A,b,Aeq,beq,uvlb,uvub,nonlcon,CtrlVar)


exitflag=nan;
output=[];





if isfield(CtrlVar.fminconUa,"Itmax")
    ItMax=CtrlVar.fminconUa.Itmax ;
else
    ItMax=10;
end

if isfield(CtrlVar.fminconUa,"TolNorm")
    TolNorm=CtrlVar.fminconUa.TolNorm ;
else
    TolNorm=1e-15;
end

if isfield(CtrlVar.fminconUa,"ReturnEachIterate")  && CtrlVar.fminconUa.ReturnEachIterate==true 
    xSeq=nan(ItMax+1,numel(x0)) ;  % only if detailed info required,
    rIterate=true;
else
    xSeq=[];
    rIterate=false;
end


JSeq=nan(ItMax+1);
JGrad=nan(ItMax+1);
NDC=blanks(ItMax) ;

gminSteepest=nan;   gammaCauchy = nan ; DeltaCauchy=0.01 ; 

x=x0;
dx=x0*0;

l=beq*0;
dl=beq*0;


It= 0;

if rIterate ; xSeq(It+1,:)=x0 ; end

while true

    It=It+1 ;

    if It>ItMax
        break
    end

    %% Newton Step

    x0=x ; %starting point
    [J0,dJdx0,H0]=func(x0) ;
    fprintf("\t Starting Point: J=%g \t norm(grad J)=%g \n",J0,norm(dJdx0))

    if ~isempty(Aeq)
        frhs=-dJdx0-Aeq'*l;
        grhs=beq-Aeq*x0;
    else
        frhs=-dJdx0;
        grhs=[];
    end

    CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical=issymmetric(H0) ;

    if CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical
        [dxNewton,dl]=solveKApeSymmetric(H0,Aeq,frhs,grhs,dx,dl,CtrlVar);
    else
        [dxNewton,dl]=solveKApe(H0,Aeq,frhs,grhs,dx,dl,CtrlVar);
    end

    xNewton=x0+dxNewton ;
    


    % Reduction as expected by a quadratic model
    dJQuadNewton=-( dJdx0'*dxNewton + dxNewton'*H0*dxNewton /2 );  % must use dJdx and H the values at start point

    [JNewton,dJdxNewton]=func(xNewton) ;
    
   


    % Actual Newton reduction:
    dJNewton=J0 - JNewton ;
    rNewton=dJNewton/dJQuadNewton ;
    fprintf("\t Full Newton step: \t J=%g \t norm(grad J)=%g \n",JNewton,norm(dJdxNewton))
    fprintf("\t \t \t \t \t (Actual Newton reduction)/(Quad reduction)=%g \n ",rNewton)


    if JNewton > J0  && CtrlVar.fminconUa.Backtracking

        fprintf("Going into Newton backtracing. \n")
        %% backtracking


        % J(x0 + gamma dx)
        %
        % dJ/dgamma = grad(J)' * dx
        %

        slope0 = -frhs'*dxNewton ;

        if slope0>0

            fprintf("slope0 not negative! \n")


        else


            Func = @(gamma)  func(x0 + gamma * dxNewton) ;

            CtrlVar.InfoLevelBackTrack=10000 ;   CtrlVar.doplots= 1 ;  CtrlVar.LineSearchAllowedToUseExtrapolation=true;
            [gmin,fmin,BackTrackInfo]=BackTracking(slope0,1,J0,JNewton,Func,CtrlVar);

            %dxNewton=gmin*dxNewton ;
            xNewton=x0+gmin*dxNewton ;
            [JNewton,dJdxNewton]=func(xNewton);
            dJNewton=J0 - JNewton ;
            dJQuad=-( dJdx0'*gmin*dxNewton + gmin*dxNewton'*H0*gmin*dxNewton /2 );  % must use dJdx and H the values at start point


        

            fprintf("\t Newton backtrack step: \t J=%g \t norm(grad J)=%g \t gmin=%g \n",JNewton,norm(dJdxNewton),gmin)
            fprintf("\t \t \t \t \t (Quad reduction)/(Actual Newton bactracking reduction)=%g \n ",dJQuad/dJNewton)



        end
    end





    %% Cauchy Step
    %
    % Minimize the quadratic model along the (neg) gradient direction.
    %
    %
    % arg min_g  Q(g) = dJdx'* g d + g d'*H* g d /2
    %
    % where d=-dJdx  and g is a scalar
    %
    %   dQ/dg= dJdx'* d + g d'*H* d =0
    %
    %       g = - ( dJdx' d) / (d' H d )
    %         =   (dJdx' dJdx)  / (dJdx' H dJdx)
    %
    %     Cauchy step is : - g dJdx
    %
    %     Cauchy degrement is : Q(g)
    %

    dJnorm = (dJdx0' * dJdx0) ;
    dJHdJ = (dJdx0' * H0 * dJdx0 ) ;  % To do, add a check if pos def, is dJdx0' H0 dJdx0 > 0

    if dJnorm <0
        fprintf("dJnorm < 0 \n")  ; % this must always be true, execpt for rounding errors.
    end

    if dJHdJ < 0
        fprintf ("dHdJ < 0 , not convex \n")
        % if H0 is nod semi-psd, the search directon is reversed, solution: pick steepest gradient as a search
        % direction with another step size,

        slope0 = -frhs'*frhs ;

        if isnan(gammaCauchy)
            gammaCauchy = -0.01 *J0/slope0 * norm(dJdx0);  % initial step size
            DeltaCauchy=gammaCauchy ;
        end


    else

        gammaCauchy = (dJdx0' * dJdx0) / (dJdx0' * H0 * dJdx0 ) ;  % To do, add a check if pos def, is dJdx0' H0 dJdx0 > 0
      

    end

    if isnan(gammaCauchy)  ||  gammaCauchy > DeltaCauchy
        gammaCauchy=DeltaCauchy;
    end

    % gammaCauchy=DeltaCauchy;  % fixing step based on quad reduction


    dxCauchy = - gammaCauchy * dJdx0/norm(dJdx0) ;
    xCauchy  =   x0 + dxCauchy ;

    [JCauchy,dJdxCauchy]=func(xCauchy) ;

    dJCauchy=J0 - JCauchy ;

    dJQuad=-( dJdx0'*dxCauchy + dxCauchy'*H0*dxCauchy /2 );  % must use dJdx and H the values at start point


    fprintf("\t Cauchy step: \t J=%g \t norm(grad J)=%g \n",JCauchy,norm(dJdxCauchy))
    rCauchy=dJCauchy/dJQuad ;
    fprintf("\t \t \t \t \t (Actual Cauchy reduction)/(Quad Reduction)=%g \n ",rCauchy)

    
    if rCauchy> 0.8
        % define step size in terms of changes in |x|

        DeltaCauchy=1.2*DeltaCauchy; % successfull
       
    elseif rCauchy < 0.5
        DeltaCauchy=DeltaCauchy/1.5; % unsuccessfull
    end



    %% Add in a line search option along steepest gradient

    Func = @(gamma)  func(x0 - gamma * dJdx0) ;
    slope0 = -frhs'*frhs ;

    % fractional reduction = -slope * step/J0 = 0.01
    % -> step = -0.01 *J0/slope

    if isnan(gminSteepest)
        b = -0.01 *J0/slope0 ;  % initial step size
    else
        b=gminSteepest;
    end

    fb=Func(b) ;
    CtrlVar.InfoLevelBackTrack=1 ;   CtrlVar.doplots= 1 ;  
    CtrlVar.LineSearchAllowedToUseExtrapolation=true;
    CtrlVar.NewtonAcceptRatio=0.01;
    [gminSteepest,fmin,BackTrackInfo]=BackTracking(slope0,b,J0,fb,Func,CtrlVar);

    xSteepest=x0-gminSteepest*dJdx0 ;
    [JSteepest,dJdxSteepest]=func(xSteepest);

    %% Compare func and quad approximation along Newton direction

    CtrlVar.fminconUa.InfoLevel=1 ;
    if CtrlVar.fminconUa.InfoLevel>=1000
        nPoints=20;
        Quad=nan(nPoints,1); step=nan(nPoints,1);   F=nan(nPoints,1);

        for II=1:nPoints
            gamma=2*II/nPoints ;
            step(II)=gamma ;

            Quad(II)=J0+( dJdx0'*gamma*dxNewton + gamma*dxNewton'*H0*gamma*dxNewton /2 );  % must use dJdx and H the values at start point
            F(II)=func(x0+gamma*dxNewton) ;
        end
        

        figFQ=FindOrCreateFigure("f(x) and f(x0)+Q(x)") ; clf(figFQ) ; 
       
       
        plot(step,Quad,"ob-",MarkerFaceColor="b")
        hold on
        plot(step,F,"*r-",MarkerFaceColor="r") 
        legend("$f(x_0)+Q(\Delta x)$","$f(x)$",interpreter="latex")

    end

    %%


    if CtrlVar.fminconUa.Step=="Auto"

        fprintf("\n J0=%g \t JNewton=%g \t JCauchy=%g \t JSteepest=%g \n \n",J0,JNewton,JCauchy,JSteepest)
        if JNewton < J0 && JNewton < JCauchy  && JNewton< JSteepest

            fprintf("Accepting Newton step \n")
            x=xNewton; J=JNewton ;  dJdx=dJdxNewton ; NDC(It)='N';

        elseif JCauchy < J0 && JCauchy < JNewton  && JCauchy < JSteepest

            fprintf("Accepting Cauchy step \n")
            x=xCauchy ;  J=JCauchy ; dJdx=dJdxCauchy ;  NDC(It)='C';

        else

            fprintf("Accepting steepest descend step \n")
            x=xSteepest ;  J=JSteepest ; dJdx=dJdxSteepest ;  NDC(It)='D';

        end


    elseif  CtrlVar.fminconUa.Step=="Newton"

        fprintf("Accepting Newton step \n")
        x=xNewton; J=JNewton ;   dJdx=dJdxNewton ;

    elseif CtrlVar.fminconUa.Step=="Cauchy"

        fprintf("Accepting Cauchy step \n")
        x=xCauchy ;J=JCauchy ; dJdx=dJdxCauchy ;

    elseif CtrlVar.fminconUa.Step=="Steepest"

        fprintf("Accepting steepest step \n")
        x=xSteepest ;J=JSteepest ; dJdx=dJdxSteepest ;

    else

        error("case not found")

    end

    if rIterate ; xSeq(It+1,:)=x ; end

    JSeq(It+1)=J;
    JGrad(It+1)=norm(dJdx); 

    if norm(dJdx) < TolNorm
        break
    end



end

if rIterate 
    output.xSeq=xSeq ;
else 
    output.xSeq=[]  ;
end

output.JSeq=JSeq;
output.JGrad=JGrad ;
output.NDC=NDC ; 

end