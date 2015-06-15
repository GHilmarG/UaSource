function fxy=FEeval(CtrlVar,connectivity,coordinates,f,xy)

ndim=2;

%
% X=a11 x + a12 y + x0
% Y=a21 x + a22 y + y0
%
% 0=a11 x(1) + a12 y(1) + x0 
% 0=a21 x(1) + a22 y(1) + y0
% 1=a11 x(2) + a12 y(2) + x0
% 0=a21 x(2) + a22 y(2) + y0
% 0=a11 x(3) + a12 y(3) + x0
% 1=a21 x(3) + a22 y(3) + y0
%
% solve for a11, a12, a21, a22 and x0 ad y0
% 6 equations with 6 unknowns
% 
%
% 0=a11 (x(1)-x0) + a12 (y(1)-y0)  
% 0=a21 (x(1)-x0) + a22 (y(1)-y0) 
% 1=a11 (x(2)-x0) + a12 (y(2)-y0) 
% 0=a21 (x(2)-x0) + a22 (y(2)-y0) 
% 0=a11 (x(3)-x0) + a12 (y(3)-y0) 
% 1=a21 (x(3)-x0) + a22 (y(3)-y0) 
%
% 6 equations with 6 unknowns
% x0=x(1) , y0=y(1);
% 4 equations with 4 unknowns
% solve using nodal coordinates
% insert into lin. transform to determin (X,Y) 
% and fxy=f_i N_i(X,Y)
%
% solution:
%  X= ( (y3-y1)(x-x1)+(x1-x3)(y-y1) )/(2A)
%  Y= ( (y1-y2)(x-x1)+(x2-x1)(y-y1) )/(2A)
%
% Hugh hassle to extend this to higher-order elements
% 

%%
clc
f=[0 0 1]'; 
coordinates=[0 0 ; 1 0 ; 0 1];
connectivity=[2 1 3];
A=TriAreaFE(coordinates,connectivity); 
Nele=1 ; nod=3; ndim=2;
fnod=reshape(f(connectivity,1),Nele,nod);
xnod=coordinates(connectivity(1,:),1);
ynod=coordinates(connectivity(1,:),2);
x=0; y=1/2;

xn=[xnod(2) xnod(1) xnod(3)];
yn=[ynod(2) ynod(1) ynod(3)];
X=( (yn(3)-yn(1))*(x-xn(1))+(xn(1)-xn(3))*(y-yn(1)) )/(2*A);
Y=( (yn(1)-yn(2))*(x-xn(1))+(xn(2)-xn(1))*(y-yn(1)) )/(2*A);
[X Y]

fun=shape_fun(1,ndim,nod,[X Y]);  
fxy=fnod*fun
%%
end

