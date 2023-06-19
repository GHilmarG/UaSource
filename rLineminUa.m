function [gammamin,rmin,BackTrackInfo] = rLineminUa(CtrlVar,UserVar,func,K,L,du0,dv0,dh0,dl0,dJdu,dJdv,dJdh,dJdl,Normalisation,M)


%  [ K  L ]  [dx0]  = [ dJdx ]
%  [ L' 0 ]  [dl0]    [ dJdl ]
%
%


% func=@(gamma,Du,Dv,Dl) CalcCostFunctionNR(UserVar,RunInfo,CtrlVar,MUA,gamma,F,fext0,L,l,cuv,Du,Dv,Dl) ;


% Newton

rNewtonFunc=@(gamma) func(gamma,du0,dv0,dl0) ;

r0=rNewtonFunc(0);
r1=rNewtonFunc(1);
slope0=-2*r0 ;



CtrlVar.uvMinimisationQuantity="Force Residuals" ;  CtrlVar.BacktracFigName="Line Search in Newton Direction" ;
[gammaminNewton,rminNewton,BackTrackInfo]=BackTracking(slope0,1,r0,r1,rNewtonFunc,CtrlVar);


%% Descent/Newton, ie using the Mass matrix as Hessian


R=[dJdu;dJdv;dJdl];
[nL,mL]=size(L);
L0=sparse(nL,nL);
% L1=speye(nL,nL); 
[nM,mM]=size(M); 
% M0=sparse(nM,nM); 
% LM0=sparse(nL,nM); 
H=[K L' ;L L0] ;
% MM=[ M M0 LM0' ; M0 M LM0' ; LM0 LM0 L1]; 


% Du=dx(1:nM) ; Dv=dx(nM+1:2*nM) ; Dl=dx(2*nM+1:end); 

Du=M\dJdu ; Dv=M\dJdv; Dl=dl0 ; % changing the search vector in uv space
rMassFunc=@(gamma) func(gamma,Du,Dv,Dl) ;

% rD0=rMassFunc(0);
% b=1 ; rDb=rMassFunc(b);
% slopeD0=-2*R'*H*(MM\R); 
% slopeD0=slope0 ;  % not quite correct 

% b = -0.1 *rD0/slopeD0 ;  % initial step size
% gamma=b ; rDb=rMassFunc(gamma);


CtrlVar.uvMinimisationQuantity="Force Residuals" ;  CtrlVar.BacktracFigName="Line Search in Mass Direction" ;
CtrlVar.LineSearchAllowedToUseExtrapolation=true;

b=1 ; slopeD0=nan ; rDb=nan ; 
[gammaminD,rminD,BackTrackInfo]=BackTracking(slopeD0,b,r0,rDb,rMassFunc,CtrlVar);


%% descent







slope0Descent=-2*R'*H*R/Normalisation ;



rDescentFunc=@(gamma) func(gamma,dJdu,dJdv,dJdl) ;


gamma=0 ; [r0,UserVar,RunInfo,rForce0,rWork0,D20,frhs,grhs]=rDescentFunc(gamma);



b = -0.1 *r0/slope0Descent ;  % initial step size
gamma=b ; rb=rDescentFunc(gamma);

CtrlVar.InfoLevelBackTrack=1000; CtrlVar.BacktracFigName="Steepest Descent" ;
CtrlVar.BacktrackingGammaMin=1e-13; CtrlVar.LineSearchAllowedToUseExtrapolation=true;
[gammaminDescent,rminDescent,BackTrackInfo]=BackTracking(slope0Descent,b,r0,rb,rDescentFunc,CtrlVar);

dNewton=[du0;dv0;dl0];
dSteepest=[dJdu;dJdv;dJdl];
p=(dNewton'*dSteepest)/(norm(dNewton)*norm(dSteepest)) ;

fprintf("angel between Newton and Steepest descent directions = %g (deg) \n",acosd(p))


%% Cauchy


R=[dJdu; dJdv ; dJdl ] ;
dx=[du0 ; dv0 ; dl0 ]  ;   % Newton direction

% Quad approximation in Newton direction
Q=@(g) (R-g*(H*dx))'*(R-g*(H*dx))/Normalisation ;

gVector=linspace(0,1.2,50)' ;
QVector=gVector*0+nan;
rVector=gVector*0+nan;

for I=1:numel(gVector)

    QVector(I)=Q(gVector(I));
    rVector(I)=rNewtonFunc(gVector(I)) ;
end

fig=FindOrCreateFigure("Q and r (Newton)") ; clf(fig) ;

hold on
plot(gVector,rVector,"o-b")
plot(gVector,QVector,"*-r")


legend("r^2","Quad approximation in Newton direction")

% Steepest Descent direction

% Quad approximation in the direction of steepest descent
dx=R ; Q=@(g) (R-g*(H*dx))'*(R-g*(H*dx))/Normalisation ;



gammaMinCauchy=R'*H*R/(R'*H'*H*R) ;

rminCauchy=rDescentFunc(gammaMinCauchy) ; 

gVector=linspace(0,2*gammaMinCauchy,50)' ;
gVector(2)=gammaMinCauchy/1000;

QVector=gVector*0+nan;
rVector=gVector*0+nan;
for I=1:numel(gVector)

    QVector(I)=Q(gVector(I));
    rVector(I)=rDescentFunc(gVector(I)) ;
end

fig=FindOrCreateFigure("Q and r steepest descent") ; clf(fig) ;

hold on
plot(gVector,rVector,"o-b")
plot(gVector,QVector,"*-r")
legend("r^2","Quad approximation in steepest descent direction")

%% Summary

fprintf("r0=%13.7g \t r1=%13.7g \t rNewton=%13.7g \t rDescent=%13.7g \t rCauchy=%13.7g \n",r0,r1,rminNewton,rminDescent,rminCauchy)
fprintf("g0=%13.8g \t g1=%13.7g \t gNewton=%13.7g \t gDescent=%13.7g \t gCauchy=%13.7g \n",0,1,gammaminNewton,gammaminDescent,gammaMinCauchy)


%%
gammamin=gammaminNewton ; rmin=rminNewton; 

return

end



















