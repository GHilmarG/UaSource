 
%% Testing lsqUa for a simple example



%  R=(x1,x2)
problemtype="[x1,x2]" ;
problemtype="[x1+x2,x2]";
problemtype="[x1^2+x2,x2]";
% problemtype="[x1^2,x2]";

isConstraint=true; 

CtrlVar.lsqUa.ItMax=20 ;
CtrlVar.lsqUa.tol=1e-30 ;
CtrlVar.lsqUa.isLSQ=true ;
CtrlVar.lsqUa.LevenbergMarquardt=1e1 ;
CtrlVar.lsqUa.Normalize=false;
CtrlVar.lsqUa.SaveIterate=true;





fun = @(x) fRK(x,problemtype)  ;
x0=[-5 ; -5] ; 

lambda= []  ;

if isConstraint
    L=[1 1 ];  c= 5  ;
else
    L=[]; c=[];
end

[xSol,lambda,resnorm,residual,g,h,output] = lsqUa(CtrlVar,fun,x0,lambda,L,c) ;




x1Vector=linspace(-10,10);
x2Vector=linspace(-10,10);
r2=nan(numel(x1Vector),numel(x2Vector));

for I=1:numel(x1Vector)
    for J=1:numel(x2Vector)
        x=[x1Vector(I) x2Vector(J)];
        R=fRK(x,problemtype);
        r2(I,J)=R'*R ;
    end
end


flsqUa=FindOrCreateFigure("lsqUa test") ; clf(flsqUa) ;
contourf(x1Vector,x2Vector,r2',10) ; axis equal tight; colorbar ; axis([-10 10 -10 10])
hold on  ;
if isConstraint
    plot(x1Vector,c-x1Vector,'r')
end

plot(xSol(1),xSol(2),'o',MarkerFaceColor='r',MarkerEdgeColor="w",MarkerSize=12)

for I=1:output.nIt

    % plot(output.xVector(1,I),output.xVector(2,I),'+r')
    text(output.xVector(1,I),output.xVector(2,I),num2str(I-1),color="r")
    
   % txt = input("RET to continue\n") ; 
    
end


flsqUaProg=FindOrCreateFigure("lsqUa progress") ; clf(flsqUaProg) ;

itVector=linspace(0,numel(output.rVector)-1);
semilogy(itVector, output.rVector,'o-')
xlabel("iteration",Interpreter="latex")
ylabel("$r^2$",Interpreter="latex")


%%
Klear

load TestUaLSQ.mat

CtrlVar.dt=0.01;



%% Standard approach
CtrlVar.InfoLevelNonLinIt=2 ; CtrlVar.doplots=1 ;
[UserVar,RunInfo,F1,l1,BCs1]=SSTREAM_TransientImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);


x1=[F1.ub;F1.vb;F1.h]; [L,cuvh,luvh1]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs1,l1);

R1=uvhRK(x1,UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
g1=-(R1 + L'*luvh1) ;
h1=-(L*x1-cuvh);
r1=full([g1;h1]'*[g1;h1]) ;



CtrlVar.uvhMatrixAssembly.ZeroFields=true;
CtrlVar.uvhMatrixAssembly.Ronly=true;
[UserVar,RunInfo,Fext0,~]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1);
Normalisation=full(Fext0'*Fext0) ;


fprintf("\t r1=%g \t Normalisation=%g \t r1/N=%g \n",r1,Normalisation,r1/Normalisation)


%
%% LSQ


fun = @(x) uvhRK(x,UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;

[L,cuvh,luvh0]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs1,l1);
if isempty(luvh0)
    luvh0=cuvh*0 ;
end

CtrlVar.lsqUa.ItMax=5 ; CtrlVar.lsqUa.tol=1e-50 ;  CtrlVar.lsqUa.isLSQ=true ;  CtrlVar.lsqUa.LevenbergMarquardt=1e30;

x0=[F1.ub;F1.vb;F1.h];
[x1,luvh1,resnorm,residual,g1,h1] = lsqUa(CtrlVar,fun,x0,luvh0,L,cuvh) ;

%% Matlab LSQ

fun = @(x) uvhRK(x,UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;

[L,cuvh,luvh0]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs1,l1);
if isempty(luvh0)
    luvh0=cuvh*0 ;
end


x0=[F1.ub;F1.vb;F1.h];
lb=[] ; ub=[] ; A=[] ; b=[] ; Aeq=L ; beq=cuvh ; 
options = optimoptions('lsqnonlin','Display','iter','MaxIterations',3);
[x1,luvh1,resnorm,residual,g1,h1] = lsqnonlin(fun,x0,lb,ub,A,b,Aeq,beq);



%%
R1=uvhRK(x1,UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
g1=-(R1 + L'*luvh1) ;
h1=-(L*x1-cuvh);
r1=full([g1;h1]'*[g1;h1]) ;
fprintf("\t r1=%g \t Normalisation=%g \t r1/N=%g \n",r1,Normalisation,r1/Normalisation)

%%

isLSQ=true ; 



[L,cuvh,luvh0]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs1,l1);
if isempty(luvh0)
    luvh0=cuvh*0 ; 
end


x0=[F1.ub;F1.vb;F1.h];
[R0,K0]=uvhRK(x0,UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;


gamma=0; fext0=1; 
dub=zeros(MUA.Nnodes,1) ; dvb=dub ; dh=dub ; dl=luvh0*0; 
[r0,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh0,cuvh,gamma,fext0) ; 

%
%   [H L' ]  [dx]  = [g]
%   [L  0 ]  [l]     [h]
%
%

H0=2*(K0'*K0) ; % 


if isLSQ
    g0=- (2*K0'*R0 + L'*luvh0 );
    h0=-(L*x0-cuvh);
else

    H0=K0 ; %
    g0=-(R0 + L'*luvh0) ;
    h0=-(L*x0-cuvh);
end


[dx,dluvh]=solveKApe(H0,L,g0,h0,x0,luvh0,CtrlVar);


x1=x0+dx ; 
n=MUA.Nnodes;
dub=dx(1:n) ;
dvb=dx(n+1:2*n) ;
dh=dx(2*n+1:3*n) ;
dl=dluvh ; 

gamma=1 ;
[r1,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh0,cuvh,gamma,fext0) ;

[R1,K1]=uvhRK(x1,UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;

luvh1=luvh0+dluvh ;

if isLSQ
    g1=- (2*K1'*R1 + L'*luvh1 );
    h1=-(L*x1-cuvh);
else
    g1=-(R1 + L'*luvh1) ;
    h1=-(L*x1-cuvh);
end

r1uvh=full(R1'*R1) ; 
r0uvh=full(R0'*R0);

r0uvhl=full([g0;h0]'*[g0;h0]) ;
r1uvhl=full([g1;h1]'*[g1;h1]) ;


fprintf("\n \n")
fprintf("\t  r0uvh=%g \t  r1uvh=%g \t   r1uvh/r0uvh=%g \n",r0uvh,r1uvh,r1uvh/r0uvh)
fprintf("\t r0uvhl=%g \t r1uvhl=%g \t r1uvhl/r0uvhl=%g \n",r0uvhl,r1uvhl,r1uvhl/r0uvhl)
fprintf("\t     r0=%g \t     r1=%g \t         r1/r0=%g \n",r0,r1,r1/r0)


fun = @(x) uvhRK(x,UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;

lambda0=luvh1;
[x1,lambda1,resnorm,residual,g1,h1] = lsqUa(CtrlVar,fun,x0,lambda0,L,cuvh) ; 







