function [x, y]=TestingGmres(CtrlVar,A,B,f,g,x0,y0)

x=[] ; y=[] ; 
tol = 1e-2;
maxit = 50;
nRestart=[];
doOnGPU=1;
%%


%%  
tgmres1=tic;
[x,flx,~,~,rvx] = gmres(A,f,nRestart,tol,maxit);
semilogy(rvx)
title('Residual Norm at Each Iteration')
tgmres1=toc(tgmres1)  ; 


%% now use equilibrate 
teq=tic;
[P,R,C] = equilibrate(A);
teq=toc(teq) ; 

if doOnGPU
    P=gpuArray(P);
    R=gpuArray(R);
    C=gpuArray(C);
end

B = R*P*A*C;
d = R*P*f;

%B=distributed(B) ; d=distributed(d) ; 

tgmres2=tic;
[y,fly,~,~,rvy] = gmres(B,d,nRestart,tol,maxit);
tgmres2=toc(tgmres2)  ; 

hold on
semilogy(rvy)
legend('Original', 'Equilibrated', 'Location', 'southeast')
title('Relative Residual Norms (No Preconditioner)')

return
%%
tgmres3=tic;
[L1,U1] = ilu(B,struct('type','ilutp','droptol',1e-2,'thresh',0));
[yp1,flyp1,~,~,rvyp1] = gmres(B,d,nRestart,tol,maxit,L1,U1);
tgmres3=toc(tgmres3)  ; 
semilogy(rvyp1)

tgmres4=tic;
[L1,U1] = ilu(B,struct('type','ilutp','droptol',1e-2,'thresh',0));
B=distributed(B) ; d=distributed(d) ;
[yp1,flyp1,~,~,rvyp1] = gmres(B,d,nRestart,tol,maxit,L1,U1);
tgmres4=toc(tgmres4)  ; 
semilogy(rvyp1)


hold off
legend('not precon', 'No precon but equilibrate','ILUTP(1e-2)','ILUTP(1e-2) dist')


fprintf(' equilibrate %5.2f sec \n',teq)
fprintf(' gmres %5.2f sec \n',tgmres1)
fprintf(' gmres %5.2f sec \n',tgmres2)
fprintf(' gmres %5.2f sec \n',tgmres3)
fprintf(' gmres %5.2f sec \n',tgmres4)
end

