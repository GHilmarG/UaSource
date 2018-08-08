
%%
clc
load TestSaveuvhAssembly

tic

for Iint=1:MUA.nip
    
    
    [Tx1,Fx1,Ty1,Fy1,Th1,Fh1,Kxu1,Kxv1,Kyu1,Kyv1,Kxh1,Kyh1,Khu1,Khv1,Khh1]=...
        uvhAssemblyIntPointImplicitSUPG(Iint,ndim,MUA,...
        bnod,hnod,unod,vnod,AGlennod,nnod,Cnod,mnod,h0nod,u0nod,v0nod,as0nod,ab0nod,as1nod,ab1nod,dadhnod,Bnod,Snod,rhonod,...
        uonod,vonod,Conod,monod,uanod,vanod,Canod,manod,...
        CtrlVar,rhow,g,Ronly,ca,sa,dt,...
        Tx0,Fx0,Ty0,Fy0,Th0,Fh0,Kxu0,Kxv0,Kyu0,Kyv0,Kxh0,Kyh0,Khu0,Khv0,Khh0);
    
    Tx=Tx+Tx1;  Fx=Fx+Fx1;
    Ty=Ty+Ty1;  Fy=Fy+Fy1;
    Th=Th+Th1;  Fh=Fh+Fh1;
    
    Kxu=Kxu+Kxu1;        Kxv=Kxv+Kxv1;
    Kyu=Kyu+Kyu1;        Kyv=Kyv+Kyv1;
    Kxh=Kxh+Kxh1;        Kyh=Kyh+Kyh1;
    Khu=Khu+Khu1;        Khv=Khv+Khv1;        Khh=Khh+Khh1;
    
end
toc


tic

parfor Iint=1:MUA.nip
    
    
    [Tx1,Fx1,Ty1,Fy1,Th1,Fh1,Kxu1,Kxv1,Kyu1,Kyv1,Kxh1,Kyh1,Khu1,Khv1,Khh1]=...
        uvhAssemblyIntPointImplicitSUPG(Iint,ndim,MUA,...
        bnod,hnod,unod,vnod,AGlennod,nnod,Cnod,mnod,h0nod,u0nod,v0nod,as0nod,ab0nod,as1nod,ab1nod,dadhnod,Bnod,Snod,rhonod,...
        uonod,vonod,Conod,monod,uanod,vanod,Canod,manod,...
        CtrlVar,rhow,g,Ronly,ca,sa,dt,...
        Tx0,Fx0,Ty0,Fy0,Th0,Fh0,Kxu0,Kxv0,Kyu0,Kyv0,Kxh0,Kyh0,Khu0,Khv0,Khh0);
    
    Tx=Tx+Tx1;  Fx=Fx+Fx1;
    Ty=Ty+Ty1;  Fy=Fy+Fy1;
    Th=Th+Th1;  Fh=Fh+Fh1;
    
    Kxu=Kxu+Kxu1;        Kxv=Kxv+Kxv1;
    Kyu=Kyu+Kyu1;        Kyv=Kyv+Kyv1;
    Kxh=Kxh+Kxh1;        Kyh=Kyh+Kyh1;
    Khu=Khu+Khu1;        Khv=Khv+Khv1;        Khh=Khh+Khh1;
    
end
toc


%%
N=10000;
A=rand(N,N);
B=rand(N,N);

%%
clc
tic
AA=distributed(A);
BB=distributed(B);
toc

tic
C=TestAddM(A,B);
t=toc;
trace(C)

tic

    CC=TestAddM(AA,BB);

tspmd=toc;
trace(CC)

fprintf('t=%f \t tspmd=%f \n',t,tspmd)













