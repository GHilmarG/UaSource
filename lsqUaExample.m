


function lsqUaExample



x0=[-10 ; 15] ;

%  R=(x1,x2)

%                                                               lsq                      H            lsq                      H
%                                                                       constraint                           unconstrained
problemtype="[x1,x2]" ;                     %                   24.5                   24.5
problemtype="[x1+x2,x2]";                   %                   25.0                    50              0                      0
problemtype="[x1^2+x2,x2]";               %                   40.915              49.999              0                      0
% problemtype="[x1^2,x2]";                  %                   16.5015             20.5917             0                      0
% problemtype="[x1^2+x2,x2^2+x1]";            %                   153.125             153.125             0                      0
% problemtype="[x1^3-100 x2,-x2^2+10 x1]" ; %                     1737.89             4052.71             0                   not conv
% problemtype="Rosenbrock" ;                  %                   1.78794              5.4718
% problemtype="lsqRosenbrock" ;      x0=[-5; -8] ;
% problemtype="[x1^2,x2^2]" ;
  % problemtype="[x1^-100 x1,0]" ;
% problemtype="[x1^-100 x1,x2^2]" ;   x0=[-5; 8] ;
problemtype="Beale" ; x0=[2 ; 0] ; 

isConstraint=false;


CtrlVar.lsqUa.ItMax=50 ;

CtrlVar.lsqUa.gTol=1e-20 ;
CtrlVar.lsqUa.dR2Tol=1e-10 ;
CtrlVar.lsqUa.dxTol=1e-20 ;


% Consistent are both:
%
% isLSQ=true and CostMeasure="r2" AND % isLSQ=false and CostMeasure="R2"
%

CtrlVar.lsqUa.isLSQ=true ; CtrlVar.lsqUa.CostMeasure="R2" ;
CtrlVar.lsqUa.isLSQ=false ; CtrlVar.lsqUa.CostMeasure="r2" ;

CtrlVar.lsqUa.LevenbergMarquardt="auto" ; % "fixed"
CtrlVar.lsqUa.LMlambda0=0 ;
CtrlVar.lsqUa.LMlambdaUpdateMethod=1 ;
CtrlVar.lsqUa.Normalize=false;
CtrlVar.lsqUa.ScaleProblem=true;
CtrlVar.lsqUa.SaveIterate=true;
CtrlVar.InfoLevelNonLinIt=1;
CtrlVar.lsqUa.Algorithm="DogLeg" ;
CtrlVar.lsqUa.DogLeg="-Newton-Cauchy-" ; 

CtrlVar.InfoLevelBackTrack=1000; 
CompareWithMatlabOpt=false;

% xSol =
%
%     1.7909
%     3.2091
%
%
% lambda =
%
%    -0.3452
%
%



fun = @(x) fRK(x,problemtype)  ;


lambda= []  ;

if isConstraint

    % x2= c - a x1
    a=1 ;
    L=[a 1 ];  c= 5  ;
    factor=1;
    L=factor*L ; c=factor*c ;

else
    L=[]; c=[];
end

[xSol,lambda,R2,r2,Slope0,dxNorm,dlambdaNorm,residual,g,h,output] = lsqUa(CtrlVar,fun,x0,lambda,L,c) ;

xSol
lambda

if numel(xSol)==2


    if problemtype=="Beale"

        xmin=-4.5; xmax=4.5 ; ymin=-4.5 ; ymax=4.5 ;
    else

        xmin= min([output.xVector(1,:) , -10])   ; ymin= min([output.xVector(2,:) , -10]) ;
        xmax= max([output.xVector(1,:) ,  10])    ; ymax= max([output.xVector(2,:) , 10]) ;

    end



    x1Vector=linspace(xmin,xmax);
    x2Vector=linspace(ymin,ymax);


    RR=nan(numel(x1Vector),numel(x2Vector));

    for I=1:numel(x1Vector)
        for J=1:numel(x2Vector)
            x=[x1Vector(I) x2Vector(J)];
            R=fRK(x,problemtype);
            RR(I,J)=R'*R ;
        end
    end


    flsqUa=FindOrCreateFigure("lsqUa test") ; clf(flsqUa) ;
    f=log10(RR'); 
    contourf(x1Vector,x2Vector,f,50,LineStyle="none") ; 
    axis equal tight; 
    cbar=colorbar ;
    title(cbar,"$\log_{10} \|R^2\|$",interpreter="latex")
    axis([xmin xmax ymin ymax])
    hold on  ;
    if isConstraint
        plot(x1Vector,c-a*x1Vector,'r')
    end

    plot(xSol(1),xSol(2),'o',MarkerFaceColor='r',MarkerEdgeColor="w",MarkerSize=12)
    plot(output.xVector(1,:),output.xVector(2,:),color="r")
    for I=1:output.nIt

        % plot(output.xVector(1,I),output.xVector(2,I),'+r')
        text(output.xVector(1,I),output.xVector(2,I),num2str(I-1),color="r")

        % txt = input("RET to continue\n") ;

    end

    
end

[flsqUaProg1,FigFound]=FindOrCreateFigure("lsqUa progress: |R| and slope") ;


ls="-"; ms="o";
if FigFound
    hold on    % to be able to compare with previous resutls, overplots and change marker
    ls="--" ;
    ms="v";
end


npoints=numel(output.R2Array);
itVector=0:npoints-1;

yyaxis left
semilogy(itVector, output.R2Array,'-',color='b',Marker=ms)
ylabel("$\|R\|^2$",Interpreter="latex")

yyaxis right
semilogy(itVector, output.r2Array,LineStyle=ls,color='r',Marker=ms)
ylabel("$\|r\|^2$",Interpreter="latex")


xlabel("iteration",Interpreter="latex")
title(sprintf("$\\|R\\|^2$ =%g, $\\|r\\|^2$=%g",R2,r2),Interpreter="latex")

[flsqUaProg,FigFound2]=FindOrCreateFigure("lsqUa progress:dx") ;

yyaxis left
semilogy(itVector+1, output.dxArray,LineStyle=ls,Color='b',Marker=ms)
ylabel("$\|\Delta x\|^2$",Interpreter="latex")


yyaxis right
semilogy(itVector,-output.Slope0Array,LineStyle=ls,Color='r',Marker=ms)
ylabel("slope",Interpreter="latex")
xlabel("iteration",Interpreter="latex")
title(sprintf("$\\|dx\\|^2$ =%g, Slope=%g",dxNorm,Slope0),Interpreter="latex")

%%






if CompareWithMatlabOpt


    % problemtype="[x1^3-100 x2,-x2^2+10 x1]" ;     L=[1 1 ];  c= 5  ; x0=[-10 ; 15] ;

    fun = @(x) fRK(x,problemtype)  ;

    lb=[] ; ub=[] ; A=[] ; b=[] ;   nonlcon=[] ;



    options = optimoptions('lsqnonlin','Display','iter','MaxIterations',30,'SpecifyObjectiveGradient',true,...
        'FunctionTolerance',1e-10,'Algorithm','interior-point',PlotFcn=@optimplotresnormUa);
    options.OptimalityTolerance = 1.000000e-20 ; options.StepTolerance = 1.000000e-20 ;
    [x1,resnorm,residual,exitflag,outputM] = lsqnonlin(fun,x0,lb,ub,A,b,L,c,nonlcon,options);

x1

end



end
