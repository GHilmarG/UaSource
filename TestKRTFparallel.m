load PlotFile

tic
[K,R,T,F]=KRTF(s,S,B,h,u,v,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,CtrlVar);
toc


Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);

FEmesh.Coordinates=coordinates;
FEmesh.Connectivity=connectivity;
FEmesh.Boundary=Boundary;
FEmesh.Nnodes=Nnodes;
FEmesh.Nele=Nele;
FEmesh.Nod=nod;
FEmesh.Nip=nip;

tic
parfor Iint=1:nip;
    [Ktest{Iint},R,T,F]=KRTFparInt(Iint,s,S,B,h,u,v,AGlen,n,C,m,FEmesh,alpha,rho,rhow,g,CtrlVar);
end
toc

