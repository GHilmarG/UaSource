





function   Fig=PlotCostVersusStepSizeAlongNewtonDirection(func,p0,dp,g0,gammaNewton,JNewton,g0SD,gammaSD,JSD)

%   PlotCostVersusStepSizeAlongNewtonDirection(func,p,dp,g0,gammaNewton,JNewton,gammaSD,JSD);

J0=func(p0);

nPoints=13;

if isnan(gammaNewton)
    gammaNewton=1;
end

gammaUp=max([1.1*gammaNewton,1]);
gammaVector=linspace(-0.1,gammaUp,nPoints);

gammaVector=unique([gammaVector 0 gammaNewton]) ; 
nPoints=numel(gammaVector);
JVector=nan(nPoints,1);

parfor I=1:nPoints

    p=p0+gammaVector(I)*dp;
    JVector(I)=func(p);

end

slope0=g0'*dp;

Fig=FindOrCreateFigure("J Newton") ; clf(Fig)
plot(gammaVector,JVector,"or-")
hold on
dgamma=0.1*gammaNewton;
plot([0 dgamma],[J0 J0+dgamma*slope0],"k--",LineWidth=2)

plot(gammaNewton,JNewton,Marker="hexagram",MarkerFaceColor="b",MarkerSize=10)
xlabel("$\gamma$",Interpreter="latex")
ylabel("$J$",Interpreter="latex")
title("Cost function ($J$) along Newton (dp) direction",Interpreter="latex")

%% gradient direction
nPoints=13;

slope0=-g0'*g0SD;
gammaUp=-0.1*J0/slope0;

if isnan(gammaSD)
    gammaSD=gammaUp;
else
    gammaUp=1.1*gammaSD;
end


gammaVector=linspace(0,gammaUp,nPoints);

gammaVector=unique([gammaVector 0 gammaSD]) ; 
nPoints=numel(gammaVector);
JVector=nan(nPoints,1);


parfor I=1:nPoints

    p=p0+gammaVector(I)*(-g0SD);
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
