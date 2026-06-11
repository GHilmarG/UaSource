function [DkvDb,DrhDb]=calcDkvDbComplex(I,s,h,u,v,coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,gfint,C,m,alpha,rho,rhow,g)
    
    
    % calculates d kv/ d b using the complex number method
    % using the fact that df/dx=Im(f(x+i dx))/dx
    
   
    
    ndim=2; dof=2; neq=dof*Nnodes;
    neqx=Nnodes ; neqy=Nnodes;
    
    N=4*nod*nod*Nele; Iind=zeros(N,1) ; Jind=Iind ; Xval=Iind;
    
    [points,weights]=sample('triangle',nip,ndim);
    % get local coordinates and weights

    %beta= sqrt(C.^(-1/m).* (sqrt(u.*u+v.*v)).^(1/m-1));
    beta2= C.^(-1/m).* (sqrt(u.*u+v.*v)).^(1/m-1) ;
    
%     dbeta=1e-10;
%     beta(I)=beta(I)+sqrt(-1)*dbeta;
    
    db=1e-15;  
    b=s-h;
    b(I)=b(I)+sqrt(-1)*db; 
    h=s-b; 
    
    %x=coordinates(:,1); y=coordinates(:,2); DT = DelaunayTri(x,y); TRI=DT.Triangulation;
    %figure(500) ; trisurf(TRI,x,y,beta) ;  title(' beta ')
    
    [kv,rh]=kvrh(s,h,coordinates,connectivity,nip,etaInt,gfint,beta2,alpha,rho,rhow,g);
        
    
    DkvDb=imag(kv)/db;
    DrhDb=imag(rh)/db;
    
end



