function [xmin,status] = CubicFit(Slope,g0,g1,g2,l1,l2,InfoLevel)

% finds the minimum based on a cubic fit given three function values
% at 0, l1, and l2 and the slope at 0, i.e.
% g0=f(0)  ; g1=f(l1) ; g2=f(l2) with Slope=f'(0)


if nargin<7
    InfoLevel=10;
end

if l1 > l2 ; ltemp=l1 ; l1=l2 ; l2=ltemp ; gtemp=g1; g1=g2 ; g2=gtemp ; end

status=0;

ab=([1/l1^2 -1/l2^2 ; -l2/l1^2 l1/l2^2 ]* [g1-Slope*l1-g0 ; g2-Slope*l2-g0])/(l1-l2);
xmin=(-ab(2)+sqrt(ab(2)^2-3*ab(1)*Slope))/3/ab(1);

if ab(1) == 0
    if InfoLevel>=10
        fprintf('CubicFit: Cubic fit not sucessfull, try parabolic fit \n ')
    end
    [ xmin , status] = parabolamin(0,l1,l2,g0,g1,g2);
    
end

if ~isreal(xmin)
    status=1;
    fprintf(' CubicFit returns an imaginary number %g+i%g \n ',real(xmin),imag(xmin))
    fprintf(' Slope=%-g \n ',Slope)
    fprintf(' %g \t %g \n',0,g0)
    fprintf(' %g \t %g \n',l1,g1)
    fprintf(' %g \t %g \n',l2,g2)
    figure ; plot([0 l1 l2],[g0 g1 g2],'-o') ; hold on ; plot([0 l1],[g0 g0+Slope*l1],'g') ; hold off
    
    if g2 < g1 && g1 < g0
        xmin=1.5*l2;
    elseif g2 > g1 && g1 > g0
        xmin=l1/2;
    else
        xmin=(l1+l2)/3;
    end
end


if isnan(xmin)
    status=1;
    fprintf(' CubicFit returns NaN \n ')
    fprintf(' Slope=%-g \n ',Slope)
    fprintf(' %g \t %g \n',0,g0)
    fprintf(' %g \t %g \n',l1,g1)
    fprintf(' %g \t %g \n',l2,g2)
end

end

