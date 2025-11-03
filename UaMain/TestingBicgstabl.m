function [x, y]=TestingBicgstabl(CtrlVar,A,B,f,g,x0,y0)

%% Start by just tesing solving A x =f 

% if trying just to solve A x = f set
n=size(A,1); 
if nargin<6 ; x0=zeros(n,1) ; end

tstart=tic ; tluinc=tic;

[L,U] = ilu(A,struct('type','ilutp','droptol',1e-6));

tluinc=toc(tluinc); disp([' ilu in  ',num2str(tluinc),' sec '])

tol=1e-7 ; maxit=50;

t1=tic ;
[sol,flag,relres,iter,resvec]=bicgstabl(A,f,tol,maxit,L,U,x0); 
t2=toc(t1);
disp([' bicgstabl  ',num2str(t2),' sec '])

x=sol(1:n) ; 
tend=toc(tstart);
disp([' total solution time  ',num2str(tend),' sec '])

figure
fprintf(' flag=%-i, iter=%-g, relres=%-g \n ',flag,iter,relres)
nnn=numel(resvec);
semilogy((0:nnn-1)/2,resvec/norm(f),'-o')
xlabel('Iteration Number')
ylabel('Relative Residual')
return

%%
tstart=tic ;
% no modification of system,


[m,n]=size(B); 
AA=[A B' ;B sparse(m,m) ] ; 
AAA=[A B' ;B sparse(1:m,1:m,1)]  ; 
bb=[f;g];
if nargin<6 ; x0=zeros(n,1) ; y0=zeros(m,1) ; end
xy0=[x0;y0];



tluinc=tic;
%setup.type = 'crout'; setup.milu = 'row'; setup.droptol = 0.1;
setup.type = 'nofill'; setup.milu = 'off';

[L1,U1] = ilu(AAA,setup);
%[L1,U1] = luinc(AA,0);
tluinc=toc(tluinc);
disp([' ilu in  ',num2str(tluinc),' sec '])



tol=1e-7 ; maxit=50;


t1=tic ;
[sol,flag,relres,iter,resvec]=bicgstabl(AA,bb,tol,maxit,L1,U1,xy0); 
t2=toc(t1);
disp([' bicgstabl  ',num2str(t2),' sec '])



x=sol(1:n) ; y=sol(n+1:end);
tend=toc(tstart);
disp([' total solution time  ',num2str(tend),' sec '])

figure
fprintf(' flag=%-i, iter=%-g, relres=%-g \n ',flag,iter,relres)
nnn=numel(resvec);
semilogy((0:nnn-1)/2,resvec/norm(bb),'-o')
xlabel('Iteration Number')
ylabel('Relative Residual')

end

