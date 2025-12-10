function [a1,a3]=CubicMelt(h1,f1,h2,f2)
    
    
 %   a3*h1^3+a1*h1=f1 ;
 %   a3*h2^3+a1*h2=f2 ;
 %  
 % Example:
 %
 % find a1 and a3 in f(x)=a1 x + a2 x^3 so that:
 %  f(1)=1 ; f(10)=-100 
 %
 %  [a1,a3]=CubicMelt(1,-1,10,-100) 
 %
 
    A=[h1^3 h1 ; h2^3 h2 ] ;
    b=[f1 ; f2]; 
    
    sol=A\b;
    a3=sol(1);
    a1=sol(2);
    
    
     h=linspace(-h2,h2);
     ab=a1*h+a3*h.^3;
     figure ; plot(h,ab) ; hold on ; plot(h,-h)
    
    
    
end