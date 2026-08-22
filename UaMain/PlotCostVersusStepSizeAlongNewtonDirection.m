





function   Fig=PlotCostVersusStepSizeAlongNewtonDirection(func,p0,dp,g0,H0,gammaNewton,JNewton,g0SD,gammaSD,JSD,gammaNewtonMax,gammaSDmax,doSteepestDecent)

%   PlotCostVersusStepSizeAlongNewtonDirection(func,p,dp,g0,gammaNewton,JNewton,gammaSD,JSD);

%%
% Quadratic approximation J(gamma)=J0+g0*dp * gamma + 0.5 dp^T H dp gamma^2 
%
% $$ Q(\gamma) = J_0 + \mathbf{g}' \Delta \mathbf{p} \, \gamma+ \frac{1}{2} (\Delta \mathbf{p})^T H\, \Delta  \mathbf{p} \, \gamma^2 $$
%
% J0 + g0'*dp *gamma + 0.5 * dp' * H0 *dp  gamma^2
%%


J0=func(p0);

nPoints=13;

if isnan(gammaNewton)
    gammaNewton=1;
end

gammaUp=min([1.1*gammaNewton,gammaNewtonMax]);
gammaDown=-0.2*gammaUp;
gammaVector=linspace(gammaDown,gammaUp,nPoints);

gammaVector=unique([gammaVector 0 gammaNewton]) ;
nPoints=numel(gammaVector);
JVector=nan(nPoints,1);

parfor I=1:nPoints

    p=p0+gammaVector(I)*dp;
    JVector(I)=func(p);

end

slope0=g0'*dp;

Fig=FindOrCreateFigure("J Newton") ; clf(Fig)
plot(gammaVector,JVector,"or-",DisplayName="$J(\gamma)$")
hold on
dgamma=0.1*gammaNewton;
plot([0 dgamma],[J0 J0+dgamma*slope0],"k--",LineWidth=2,DisplayName="slope at origin")

plot(gammaNewton,JNewton,Marker="hexagram",MarkerFaceColor="b",LineStyle="none",MarkerSize=10,DisplayName="Minimum as found")
xlabel("$\gamma$",Interpreter="latex")
ylabel("$J$",Interpreter="latex")
title("Cost function ($J$) along Newton (dp) direction",Interpreter="latex")
subtitle(sprintf("$J(\\gamma)$=%g  and $\\gamma$=%g",JNewton,gammaNewton),Interpreter="latex") 

gammaRange=linspace(min(gammaVector),max(gammaVector)) ;
Q= J0 + (g0'*dp) *gammaRange + (0.5 * dp' * H0 *dp ) * gammaRange.^2 ;
hold on ; plot(gammaRange,Q,"--",color="b",DisplayName="Local quadradic approximation")

lg=legend(Interpreter="latex");
%% gradient direction


if doSteepestDecent

    nPoints=13;

    slope0=g0'*g0SD;
    gammaUp=-0.1*J0/slope0;

    if isnan(gammaSD)
        gammaSD=gammaUp;
    else
        gammaUp=1.1*gammaSD;
    end

    gammaUp=min(gammaUp,gammaSDmax) ;

    gammaVector=linspace(0,gammaUp,nPoints);

    gammaVector=unique([gammaVector 0 gammaSD]) ;
    nPoints=numel(gammaVector);
    JVector=nan(nPoints,1);


    parfor I=1:nPoints

        p=p0+gammaVector(I)*g0SD;
        JVector(I)=func(p);

    end



    Fig=FindOrCreateFigure("J grad") ; clf(Fig)
    plot(gammaVector,JVector,"or-")
    hold on
    dgamma=0.1*gammaUp;
    plot([0 dgamma],[J0 J0+dgamma*slope0],"k--",LineWidth=2)

    plot(gammaSD,JSD,Marker="hexagram",MarkerFaceColor="b",MarkerSize=10)

    xlabel("$\gamma$",Interpreter="latex")
    ylabel("$J$",Interpreter="latex")
    title("Cost function ($J$) along negative grad direction",Interpreter="latex")
    subtitle(sprintf("$J(\\gamma)$=%g  and  $\\gamma$=%g",JSD,gammaSD),Interpreter="latex")
end
%%




alpha=angleBetweenVector(-g0SD,dp);

fprintf("angle between dp and (-g0) is %f degrees\n",alpha)

%%
end


function theta = angleBetweenVector(x,y)

if all(x==0) || all(y==0)
    theta = NaN;
    return
end
a = x*norm(y) - y*norm(x);
b = x*norm(y) + y*norm(x);
theta = 2 * atan2d( norm(a),norm(b) ); % degrees

end
