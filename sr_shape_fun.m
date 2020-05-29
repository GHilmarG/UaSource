function fun = sr_shape_fun(B,nod)

[m,~] = size(B);
fun=zeros(m,nod);

 one=1.0 ; two=2.0 ; d3=3.0 ; d4=4.0 ;  d9=9.0;
 d27=27.0 ; 


c1=B(:,1);
c2 =B(:,2);
c3 = B(:,3);

switch nod
    case 3 %linear triangle
        fun = [c1 c2 c3]  ;
    case 6 %quadratic triangle
        fun(:,1)=(two.*c1-one).*c1;
        fun(:,2)=d4.*c2.*c1;
        fun(:,3)=(two.*c2-one).*c2 ;
        fun(:,4)=d4.*c2.*c3      ;
        fun(:,5)=(two.*c3-one).*c3;
        fun(:,6)=d4.*c1.*c3 ;
    case 10 %cubic triangle
        fun(:,1)= ((d3*c1-one)*(d3*c1-two)*c1)/two;
        fun(:,2)= -(d9*(d3*c1-one)*(c1+c3-one)*c1)/two;
        fun(:,3)=  (d9*(d3*c1+d3*c3-two)*(c1+c3-one)*c1)/two;
        fun(:,4)=-((d3*c1+d3*c3-one)*(d3*c1+d3*c3-two)*(c1+c3-one))/two    ;
        fun(:,5)=  (d9*(d3*c1+d3*c3-two)*(c1+c3-one)*c3)/two;
        fun(:,6)= -(d9*(c1+c3-one)*(d3*c3-one)*c3)/two;
        fun(:,7)= ((d3*c3-one)*(d3*c3-two)*c3)/two;
        fun(:,8)=  (d9*(d3*c3-one)*c1*c3)/two;
        fun(:,9)=  (d9*(d3*c1-one)*c1*c3)/two;
        fun(:,10)=-d27*((c3-one)+c1)*c1*c3;
    otherwise
        disp('this type of triangle is not compatible with mapping function')
        
end