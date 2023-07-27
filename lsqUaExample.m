


function lsqUaExample



%  R=(x1,x2)

%                                                               lsq                      H            lsq                      H
%                                                                       constraint                           unconstraint
problemtype="[x1,x2]" ;                     %                   24.5                   24.5
% problemtype="[x1+x2,x2]";                   %                   25.0                    50              0                      0
% problemtype="[x1^2+x2,x2]";               %                   40.915              49.999              0                      0
% problemtype="[x1^2,x2]";                  %                   16.5015             20.5917             0                      0
% problemtype="[x1^2+x2,x2^2+x1]";            %                   153.125             153.125             0                      0
% problemtype="[x1^3-100 x2,-x2^2+10 x1]" ; %                     1737.89             4052.71             0                   not conv
problemtype="Rosenbrock" ;                %                   1.78794              5.4718


isConstraint=false;


CtrlVar.lsqUa.ItMax=20 ;

CtrlVar.lsqUa.gTol=1e-20 ;
CtrlVar.lsqUa.dR2Tol=1e-2 ;
CtrlVar.lsqUa.dxTol=1e-20 ;

CtrlVar.lsqUa.isLSQ=true ;
CtrlVar.lsqUa.LevenbergMarquardt="auto" ; % "fixed"
CtrlVar.lsqUa.LMlambda0=0 ;
CtrlVar.lsqUa.LMlambdaUpdateMethod=1 ;
CtrlVar.lsqUa.Normalize=false;
CtrlVar.lsqUa.ScaleProblem=true;
CtrlVar.lsqUa.SaveIterate=true;

CtrlVar.lsqUa.Algorithm="DogLeg" ;

CompareWithMatlabOpt=true;


fun = @(x) fRK(x,problemtype)  ;
x0=[-10 ; 15] ;
% x0=[5; 0] ;


lambda= []  ;

if isConstraint

    % x2= c - a x1
    a=1 ; c=5 ; 
    L=[a 1 ];  c= 5  ;

else
    L=[]; c=[];
end

[xSol,lambda,R2,Slope0,g2,residual,g,h,output] = lsqUa(CtrlVar,fun,x0,lambda,L,c) ;


xmin=min(xSol(1)-1,-10) ; ymin=min(xSol(2)-1,-10) ;
xmax=max(xSol(1)+1,10) ; ymax=max(xSol(2)+1,10) ;

x1Vector=linspace(xmin,xmax);
x2Vector=linspace(ymin,ymax);


r2=nan(numel(x1Vector),numel(x2Vector));

for I=1:numel(x1Vector)
    for J=1:numel(x2Vector)
        x=[x1Vector(I) x2Vector(J)];
        R=fRK(x,problemtype);
        r2(I,J)=R'*R ;
    end
end


flsqUa=FindOrCreateFigure("lsqUa test") ; clf(flsqUa) ;
contourf(x1Vector,x2Vector,r2',20) ; axis equal tight; colorbar ; axis([xmin xmax ymin ymax])
hold on  ;
if isConstraint
    plot(x1Vector,c-a*x1Vector,'r')
end

plot(xSol(1),xSol(2),'o',MarkerFaceColor='r',MarkerEdgeColor="w",MarkerSize=12)

for I=1:output.nIt

    % plot(output.xVector(1,I),output.xVector(2,I),'+r')
    text(output.xVector(1,I),output.xVector(2,I),num2str(I-1),color="r")

    % txt = input("RET to continue\n") ;

end


[flsqUaProg,FigFound]=FindOrCreateFigure("lsqUa progress") ;


if FigFound
    hold on
end


npoints=numel(output.R2Array);
itVector=0:npoints-1;
yyaxis right
semilogy(itVector, output.g2Array,'o-')
semilogy(itVector,-output.Slope0Array,'o-')
ylabel("slope",Interpreter="latex")
yyaxis left
semilogy(itVector, output.R2Array,'o-')


ylabel("$\|R\|^2$",Interpreter="latex")
xlabel("iteration",Interpreter="latex")
title(sprintf("$\\|R\\|^2$ =%g, slope=%g, $\\|g\\|^2$=%g",R2,Slope0,g2),Interpreter="latex")

%%

if CompareWithMatlabOpt


    % problemtype="[x1^3-100 x2,-x2^2+10 x1]" ;     L=[1 1 ];  c= 5  ; x0=[-10 ; 15] ;

    fun = @(x) fRK(x,problemtype)  ;

    lb=[] ; ub=[] ; A=[] ; b=[] ;   nonlcon=[] ;



    options = optimoptions('lsqnonlin','Display','iter','MaxIterations',30,'SpecifyObjectiveGradient',true,...
        'FunctionTolerance',1e-10,'Algorithm','interior-point',PlotFcn=@optimplotresnormUa);

    [x1,resnorm,residual,exitflag,outputM] = lsqnonlin(fun,x0,lb,ub,A,b,L,c,nonlcon,options);



end



end
