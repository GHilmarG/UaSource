function [fun]=shape_fun(Iint,ndim,nod,points)
    
    %   computes the values of the shape functions.
    %   to local coordinates
    
    fun=zeros(nod,1) ;
    
    %xii=zeros(20,1) ; etai=zeros(20,1) ; zetai=zeros(20,1) ;
    pt125=0.125 ; pt25=0.25 ; pt5=0.5 ; pt75=0.75 ; one=1.0 ; two=2.0 ; d3=3.0 ; d4=4.0 ; d8=8.0 ; d9=9.0;
    d16=16.0 ; d27=27.0 ; d32=32.0 ; d64=64.0 ; d128=128.0;
    
    
    switch ndim
        case 1 % one dimensional case
            xi=points(Iint,1);
            switch nod
                case 2
                    t1=-one-xi ;
                    t2= one-xi;
                    fun(1)=t2/two;
                    fun(2)=-t1/two;
                case 3
                    t1=-one-xi ;
                    t2=-xi ;
                    t3=one-xi;
                    fun(1)=t2*t3/two ;
                    fun(2)=-t1*t3 ;
                    fun(3)=t1*t2/two;
                case 4
                    t1=-one-xi ;
                    t2=-one/d3-xi ;
                    t3=one/d3-xi ;
                    t4=one-xi;
                    fun(1)=t2*t3*t4*d9/d16  ;
                    fun(2)=-t1*t3*t4*d27/d16;
                    fun(3)=t1*t2*t4*d27/d16 ;
                    fun(4)=-t1*t2*t3*d9/d16;
                case 5
                    t1=-one -xi ;
                    t2=-pt5-xi ;
                    t3=-xi ;
                    t4=pt5-xi ;
                    t5=one-xi;
                    fun(1)=t2*t3*t4*t5*two/d3 ;
                    fun(2)=-t1*t3*t4*t5*d8/d3 ;
                    fun(3)=t1*t2*t4*t5*d4 ;
                    fun(4)=-t1*t2*t3*t5*d8/d3;
                    fun(5)=t1*t2*t3*t4*two/d3;
                otherwise
                    disp('wrong number of nodes in shape_fun')
            end
        case 2 % two dimensional case
            c1=points(Iint,1);
            c2=points(Iint,2);
            c3=one-c1-c2 ;
            xi=points(Iint,1);
            eta=points(Iint,2);
            etam=pt25*(one-eta);
            etap=pt25*(one+eta);
            xim=pt25*(one-xi);
            xip=pt25*(one+xi);
            switch nod
                case 3
                    fun = [c1 ;c3 ;c2]  ;
                case 6
                    fun(1)=(two*c1-one)*c1;
                    fun(2)=d4*c3*c1;
                    fun(3)=(two*c3-one)*c3 ;
                    fun(4)=d4*c2*c3      ;
                    fun(5)=(two*c2-one)*c2;
                    fun(6)=d4*c1*c2 ;
                case 10
                    fun(1)= ((d3*c1-one)*(d3*c1-two)*c1)/two;
                    fun(2)= -(d9*(d3*c1-one)*(c1+c2-one)*c1)/two;
                    fun(3)=  (d9*(d3*c1+d3*c2-two)*(c1+c2-one)*c1)/two;
                    fun(4)=-((d3*c1+d3*c2-one)*(d3*c1+d3*c2-two)*(c1+c2-one))/two    ;
                    fun(5)=  (d9*(d3*c1+d3*c2-two)*(c1+c2-one)*c2)/two;
                    fun(6)= -(d9*(c1+c2-one)*(d3*c2-one)*c2)/two;
                    fun(7)= ((d3*c2-one)*(d3*c2-two)*c2)/two;
                    fun(8)=  (d9*(d3*c2-one)*c1*c2)/two;
                    fun(9)=  (d9*(d3*c1-one)*c1*c2)/two;
                    fun(10)=-d27*((c2-one)+c1)*c1*c2;
                case 15
                    t1=c1-pt25  ;
                    t2=c1-pt5 ;
                    t3=c1-pt75   ;
                    t4=c2-pt25;
                    t5=c2-pt5   ;
                    t6=c2-pt75 ;
                    t7=c3-pt25  ;
                    t8=c3-pt5 ;
                    t9=c3-pt75;
                    fun(1)=d32/d3*c1*t1*t2*t3   ;
                    fun(2)=d128/d3*c3*c1*t1*t2;
                    fun(3)=d64*c3*c1*t1*t7      ;
                    fun(4)=d128/d3*c3*c1*t7*t8;
                    fun(5)=d32/d3*c3*t7*t8*t9   ;
                    fun(6)=d128/d3*c2*c3*t7*t8;
                    fun(7)=d64*c2*c3*t4*t7      ;
                    fun(8)=d128/d3*c2*c3*t4*t5;
                    fun(9)=d32/d3*c2*t4*t5*t6   ;
                    fun(10)=d128/d3*c1*c2*t4*t5;
                    fun(11)=d64*c1*c2*t1*t4     ;
                    fun(12)=d128/d3*c1*c2*t1*t2;
                    fun(13)=d128*c1*c2*t1*c3    ;
                    fun(15)=d128*c1*c2*c3*t4;
                    fun(14)=d128*c1*c2*c3*t7      ;
                case 4
                    fun=[d4*xim*etam;d4*xim*etap;d4*xip*etap;d4*xip*etam];
                case 5
                    fun=[d4*xim*etam-pt25*(one-xi.^2)*(one-eta.^2);			  ...
                        d4*xim*etap-pt25*(one-xi.^2)*(one-eta.^2);			  ...
                        d4*xip*etap-pt25*(one-xi.^2)*(one-eta.^2);			  ...
                        d4*xip*etam-pt25*(one-xi.^2)*(one-eta.^2);			  ...
                        (one-xi.^2)*(one-eta.^2)];
                case 8
                    fun=[d4*etam*xim*(-xi-eta-one);d32*etam*xim*etap;                   ...
                        d4*etap*xim*(-xi+eta-one);d32*xim*xip*etap;                    ...
                        d4*etap*xip*(xi+eta-one); d32*etap*xip*etam;                   ...
                        d4*xip*etam*(xi-eta-one); d32*xim*xip*etam];
                case 9
                    etam=eta-one;
                    etap=eta+one;
                    xim=xi-one;
                    xip=xi+one;
                    fun=[pt25*xi*xim*eta*etam;-pt5*xi*xim*etap*etam;                   ...
                        pt25*xi*xim*eta*etap;-pt5*xip*xim*eta*etap;                    ...
                        pt25*xi*xip*eta*etap;-pt5*xi*xip*etap*etam;                    ...
                        pt25*xi*xip*eta*etam;-pt5*xip*xim*eta*etam;                    ...
                        xip*xim*etap*etam];
                otherwise
                    disp('wrong number of nodes in shape_fun')
            end
        case 3 % 3 dimensional case
            xi=points(Iint,1);
            eta=points(Iint,2);
            zeta=points(Iint,3);
            etam=one-eta ;
            xim=one-xi  ;
            zetam=one-zeta;
            etap=eta+one ;
            xip=xi+one   ;
            zetap=zeta+one;
            switch nod
                case 4
                    fun(1)=xi   ;
                    fun(2)=eta ;
                    fun(3)=zeta ;
                    fun(4)=one-fun(1)-fun(2)-fun(3);
                case 8
                    fun=[pt125*xim*etam*zetam;pt125*xim*etam*zetap;                     ...
                        pt125*xip*etam*zetap;pt125*xip*etam*zetam;                     ...
                        pt125*xim*etap*zetam;pt125*xim*etap*zetap;                     ...
                        pt125*xip*etap*zetap;pt125*xip*etap*zetam];
                case 14  %  type 6 elemen
                    fun(1) = (xi*eta+xi*zeta+two*xi+eta*zeta+two*eta+two*zeta+two)*      ...
                        (xi-one)*(eta-one)*(zeta-one)/d8;
                    fun(2) =-(xi*eta-xi*zeta+two*xi-eta*zeta+two*eta-two*zeta+two)*      ...
                        (xi-one)*(eta-one)*(zeta+one)/d8;
                    fun(3) =-(xi*eta-xi*zeta+two*xi+eta*zeta-two*eta+two*zeta-two)*      ...
                        (xi+one)*(eta-one)*(zeta+one)/d8;
                    fun(4) = (xi*eta+xi*zeta+two*xi-eta*zeta-two*eta-two*zeta-two)*      ...
                        (xi+one)*(eta-one)*(zeta-one)/d8;
                    fun(5) =-(xi+one)*(xi-one)*(eta-one)*(zeta+one)*(zeta-one)/two;
                    fun(6) =-(xi-one)*(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two;
                    fun(7) = (xi+one)*(xi-one)*(eta+one)*(eta-one)*(zeta+one)/two;
                    fun(8) = (xi+one)*(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two;
                    fun(9) =-(xi+one)*(xi-one)*(eta+one)*(eta-one)*(zeta-one)/two;
                    fun(10)= (xi*eta-xi*zeta-two*xi+eta*zeta+two*eta-two*zeta-two)*      ...
                        (xi-one)*(eta+one)*(zeta-one)/d8;
                    fun(11)=-(xi*eta+xi*zeta-two*xi-eta*zeta+two*eta+two*zeta-two)*      ...
                        (xi-one)*(eta+one)*(zeta+one)/d8;
                    fun(12)=-(xi*eta+xi*zeta-two*xi+eta*zeta-two*eta-two*zeta+two)*      ...
                        (xi+one)*(eta+one)*(zeta+one)/d8;
                    fun(13)= (xi*eta-xi*zeta-two*xi-eta*zeta-two*eta+two*zeta+two)*      ...
                        (xi+one)*(eta+one)*(zeta-one)/d8;
                    fun(14)= (xi+one)*(xi-one)*(eta+one)*(zeta+one)*(zeta-one)/two;
                case 20
                    xii=[-1;-1;-1;0;1;1;1;0;-1;-1;1;1;-1;-1;-1;0;1;1;1;0];
                    etai=[-1;-1;-1;-1;-1;-1;-1;-1;0;0;0;0;1;1;1;1;1;1;1;1];
                    zetai=[-1;0;1;1;1;0;-1;-1;-1;1;1;-1;-1;0;1;1;1;0;-1;-1];
                    for l=1:20;
                        xi0=xi*xii(l);
                        eta0=eta*etai(l);
                        zeta0=zeta*zetai(l);
                        if l==4 || l==8 || l==16 ||l==20
                            fun(l)=pt25*(one-xi*xi)*(one+eta0)*(one+zeta0);
                        elseif l>=9 && l<=12
                            fun(l)=pt25*(one+xi0)*(one-eta*eta)*(one+zeta0);
                        elseif l==2 || l==6 || l==14 || l==18
                            fun(l)=pt25*(one+xi0)*(one+eta0)*(one-zeta*zeta);
                        else
                            fun(l)=pt125*(one+xi0)*(one+eta0)*(one+zeta0)*(xi0+eta0+zeta0-2);
                        end
                    end
                otherwise
                    disp('wrong number of nodes in shape_fun')
            end
        otherwise
            disp('wrong number of dimensions in shape_fun')
    end
end
