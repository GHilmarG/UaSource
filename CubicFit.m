function [xmin,status] = CubicFit(Slope,y0,y1,y2,x1,x2,InfoLevel)

% finds the minimum based on a cubic fit given three function values
% at 0, l1, and l2 and the slope at 0, i.e.
% g0=f(0)  ; g1=f(l1) ; g2=f(l2) with Slope=f'(0)


if nargin<7
    InfoLevel=0;
end


if any(isnan([y0;y1;y2;x1;x2]))

    fprintf("CubicFit: x1=%g \t x2=%g \t y0=%g \t y1=%g \t y2=%g \n",x1,x2,y0,y1,y2)
    error("CubicFit:nanOnInput","input values contain nan")

end

% If not correctly ordered on input, make sure that 0 < x1 < x2
if x1 > x2 ; ltemp=x1 ; x1=x2 ; x2=ltemp ; gtemp=y1; y1=y2 ; y2=gtemp ; end


status=0; 


% First check some special cases for which fitting a cubic polynomial can not be done.
if x1==x2

    % If on input x1=x2, find a parabolic min
    
    D=((y1-y0)/x1-Slope) ;

    if D~=0

        % This is a parabolic minimum based on slope at x=x0=0 and the value
        % pairs (0,y0) and (x1,y1)

        xmin=-x1*Slope/2/D;

    else

        % OK now, both x1 and x2 are same on input, and slope at origin (ie x=0) is
        % same as the slope from x1 to x2. Just return a value close to, but smaller than x1

        xmin=x1*0.9; 

    end

    if ~isnan(xmin)
        status=1;
    end

else

    % OK, here we are in the situation for which this routine is really written for.

    % This is the min of a cubic polynomial going through (0,y0), (x1,y1), (x2,y2) and with a slope Slope0 at x=0

    ab=([1/x1^2 -1/x2^2 ; -x2/x1^2 x1/x2^2 ]* [y1-Slope*x1-y0 ; y2-Slope*x2-y0])/(x1-x2);
    xmin=(-ab(2)+sqrt(ab(2)^2-3*ab(1)*Slope))/3/ab(1);


    % It can quite easily happen that the cubic fit does not have a minimum.
    % In which case a parabolic fit is a better option to try.

    if ab(1) == 0 || ~isreal(xmin) || isnan(xmin)

        if ~isnan(Slope0)
            D=((y1-y0)/x1-Slope) ;
            xmin=-x1*Slope/2/D;
            if ~isnan(xmin)
                status=1;
            end

        else
            % No slope information available, just fit a parabola through (0,y0), (x1,y1), and (x2,y2)
            [ xmin , status] = parabolamin(0,x1,x2,y0,y1,y2);
        end

    end

end

if isnan(xmin) || ~isreal(xmin)

    % so 0 < x1 < x2 and there are no nans in input, and I still get nan for xmin.
    % Simply use bisection

    if y0 < y1 && y1 < y2
        xmin=x1/2;

    elseif y0 > y1 && y1 > y2

        xmin=x1+(x2-x1)/2 ;

    elseif y0 < y1 && y1 > y2

        if y2 < y0
            xmin=x1+(x2-x1)/2 ;
        else
            xmin=x1/2;
        end

    end

if InfoLevel>10000


    FindOrCreateFigure("CubicFit") ;   plot([0 x1 x2],[y0 y1 y2],'o') ; hold on ; 
    xline(xmin,"--",Label="$x_{\mathrm{min}}$",Interpreter="latex") ; 
    plot([0 x1/10],[y0 y0+Slope*x1/10],'r') ; 
    xlabel("$x$",interpreter="latex") ;  
    ylabel("$y$",interpreter="latex") ;

        

end



end

