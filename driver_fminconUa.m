

doPlotFunction=1;


a=1 ; 



func=@(x) RosenbrockFunction(x) ;




%% optimisation

A=[] ; b=[] ; Aeq=[] ; beq=[] ; xLb=[] ; xUb=[] ; nonlcon=[] ;

x0=[-1 ; 6] ;
% x0=[-3 ;-1] ;


CtrlVar.fminconUa.Itmax=10;




CtrlVar.fminconUa.TolNorm=1e-6 ;
CtrlVar.fminconUa.Step="Newton" ;
% CtrlVar.fminconUa.Step="Cauchy" ;
% CtrlVar.fminconUa.Step="Steepest" ;
% CtrlVar.fminconUa.Step="Auto" ;

CtrlVar.fminconUa.ReturnEachIterate=true ; 
CtrlVar.fminconUa.Backtracking=false ;
[x,Jexit,exitflag,output] = fminconUa(func,x0,A,b,Aeq,beq,xLb,xUb,nonlcon,CtrlVar) ;



%% Plotting


if doPlotFunction
    n=100 ;
    x1=linspace(-3,2,n);
    x2=linspace(-1,7,n);

    [X,Y]=ndgrid(x1,x2); Z=X*0; L=X*0 ;   % for plotting

    for I=1:n
        for J=1:n
            x=[x1(I) , x2(J)] ;
            [f,grad,Hess]=func(x);
            Z(I,J)=f;
            lmin=eigs(Hess,1,'smallestreal');
            L(I,J)=lmin;
        end
    end


    fpath=FindOrCreateFigure("path") ; clf(fpath);
    contourf(X,Y,Z,40) ; colorbar ; axis equal
    hold on
    xSeq=output.xSeq;


    if ~isempty(xSeq)
        plot(xSeq(:,1),xSeq(:,2),'-or',MarkerFaceColor="r")
    end
    xlabel("$x$",Interpreter="latex") ; ylabel("$y$",Interpreter="latex")

    feig=FindOrCreateFigure("min eig") ; clf(feig);
    contourf(X,Y,L,40) ; colorbar ; axis equal ; ModifyColormap(0,20);
    title("Smallest eigenvalue")
    
    xlabel("$x$",Interpreter="latex") ; ylabel("$y$",Interpreter="latex")



end