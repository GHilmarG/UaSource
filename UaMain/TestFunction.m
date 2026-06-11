
function  [J,dJdp,Hess,fOuts]=TestFunction(p,TestFunctionParameters)
x=p(1);
y=p(2);
a=TestFunctionParameters.a;
b=TestFunctionParameters.b;
switch lower(TestFunctionParameters.Name)
    
    
    case 'sphere'
        
        
        J=p(1)^2+p(2)^2 ;
        dJdp=[2*p(1) ; 2*p(2)];
        
    case 'ellipse'
        
  
        
        J=(x./a).^2+(y/b).^2 ;
        dJdp=[2*x./a ; 2*y./b];
        
    case 'rosenbrock'
        
        J=b*(y-x^2)^2+(x-a)^2;
        dJdp=[-4*b*(y-x^2)*x+2*(x-a) ; 2*b*(y-x^2)];
        
    case 'matyas'
        
       J= a*(x^2+y^2)-b*x*y;
       dJdp=[2*a*x - b* y ; 2*a*y - b*x];
        
    otherwise
        
        error('sdfa')
        
end


Hess=1;
fOuts.MisfitOuts.I=NaN;
fOuts.RegOuts.R=NaN;

end