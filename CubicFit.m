function [xmin,status] = CubicFit(Slope,y0,y1,y2,x1,x2,InfoLevel)

% finds the minimum based on a cubic fit given three function values
% at 0, l1, and l2 and the slope at 0, i.e.
% g0=f(0)  ; g1=f(l1) ; g2=f(l2) with Slope=f'(0)


if nargin<7
    InfoLevel=0;
end

if x1 > x2 ; ltemp=x1 ; x1=x2 ; x2=ltemp ; gtemp=y1; y1=y2 ; y2=gtemp ; end

status=0;

ab=([1/x1^2 -1/x2^2 ; -x2/x1^2 x1/x2^2 ]* [y1-Slope*x1-y0 ; y2-Slope*x2-y0])/(x1-x2);
xmin=(-ab(2)+sqrt(ab(2)^2-3*ab(1)*Slope))/3/ab(1);



% f(x)= a0 + a1 x + a2 x^2 + a3 x^3
%
% f(0)= y0 = a0
% df/dx|_0= =Slope = a1
%
% f(l1)=l1 = a0 * a1 l1 + a2 l1^2 + a3 * l1^3
% f(l2)=l2 = a0 * a1 l2 + a2 l2^2 + a3 * l2^3
%


% It can quite easity happen that the cubic fit does not have a minimum.
% In which case a parabolic fit is a better option to try.

if ab(1) == 0 || ~isreal(xmin)
    if InfoLevel>=10
        fprintf('CubicFit: Minima not found within the range [0,%f]. Minimiser at one of the endpoints.  \n ',x2)
    end
    
    [ xmin , status] = parabolamin(0,x1,x2,y0,y1,y2);

    
end




if ~isreal(xmin) || InfoLevel>10000 

 % Cubic fit
    A=[ 1  0  0   0       ; ...
        1  x1  x1^2  x1^3 ; ...
        1  x2  x2^2  x2^3 ; ...
        0   1   0     0    ] ;

    b=[y0 ; y1 ; y2 ; Slope] ;
    sol=A\b ; a0=sol(1) ; a1=sol(2) ; a2=sol(3) ; a3=sol(4) ;

 
    % least-squares parabolic fit
 A=   [ 1   0   0    ; ...
        1  x1  x1^2  ; ...
        1  x2  x2^2  ; ...
        0   1   0    ] ;

    b=[y0 ; y1 ; y2 ; Slope] ;
    warning('off','MATLAB:rankDeficientMatrix')
    sol=A\b ; A0=sol(1) ; A1=sol(2) ; A2=sol(3) ; 
  
    fprintf(' CubicFit returns an imaginary number %g+i%g \n ',real(xmin),imag(xmin))
    fprintf(' Slope=%-g \n ',Slope)
    fprintf(' %g \t %g \n',0,y0)
    fprintf(' %g \t %g \n',x1,y1)
    fprintf(' %g \t %g \n',x2,y2)

    FindOrCreateFigure("CubicFit")
    
    plot([0 x1 x2],[y0 y1 y2],'o') ;
    hold on ; 
    plot([0 x1],[y0 y0+Slope*x1],'g') ;
    x=linspace(0,2*x2) ;
    y=a0+a1*x+a2*x.^2+a3*x.^3 ;
    plot(x,y,LineWidth=2,LineStyle="--")
    y=A0+A1*x+A2*x.^2 ;
    plot(x,y)

    if ~isnan(xmin)
        xline(xmin,'--') ;
    end
    hold off

    legend("data","slope","cubic fit","least squares parabolic fit","xmin estimated")



 
end



end

