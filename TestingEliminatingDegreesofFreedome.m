load TestSave A B f g x0 y0 CtrlVar ;

n=size(A,1) ; m=size(B,1);

tic

C=sparse(m,m);
AA=[A B' ;B -C] ; bb=[f;g];
sol=AA\bb; toc
x=sol(1:n) ; y=sol(n+1:n+m);


tic
[I,iConstrainedDOF]=ind2sub(size(B),find(B==1)); iConstrainedDOF=iConstrainedDOF(:);
iFreeDOF=setdiff(1:n,iConstrainedDOF); iFreeDOF=iFreeDOF(:);


AA=A; ff=f;
AA(iConstrainedDOF,:)=[]; AA(:,iConstrainedDOF)=[]; ff(iConstrainedDOF)=[];
sol2=AA\ff; 

xx=zeros(n,1) ; 
xx(iConstrainedDOF)=g ; 
xx(iFreeDOF)=sol2;
yy=B*(f-A*xx);
toc






