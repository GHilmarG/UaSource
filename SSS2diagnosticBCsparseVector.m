

function SSTREAM2d(s,h,u,v,coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,gfint,AGlen,C,L,Lb,lambdau,n,m,alpha,rho,rhow,g,Itime)
    
    fprintf(' calculating kv \n')
    
    %if any(h<0) ; error(' thickness negative ') ; end
    beta2= C.^(-1/m).* (sqrt(u.*u+v.*v)).^(1/m-1) ;
    
    [kv,rh]=kvrh(s,h,coordinates,connectivity,nip,etaInt,gfint,beta2,alpha,rho,rhow,g)
    