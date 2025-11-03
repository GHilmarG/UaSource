
function [ xmin , status] = parabolamin(a,b,c,fa,fb,fc,InfoLevel)

if nargin <7
    InfoLevel=10;
end

% [ xmin , status] = parabolamin(a,b,c,fa,fb,fc)
% finds the abscissa x that is the minimum of a parabola trhough three points
% status:  0 if everything went OK
%          1 problem, for example: not convex, input values contains NaN
%

status=0;  % 0 is everything went OK as far as can be seen


if nargin==2
    xvec=a ; fvec=b;
    if numel(xvec) ~=numel(fvec) ; error(' vectors of unequal length ') ; end
    if numel(xvec) ~=3           ; error(' incorrect number of elements on input') ; end
    a=xvec(1) ; b=xvec(2) ; c=xvec(3);
    fa=fvec(1) ; fb=fvec(2) ; fc=fvec(3);
end

if c < b
    ctemp=c ; fctemp=fc ; 
    c=b ; fc=fb ;
    b=ctemp ; fb=fctemp ; 
end


if ~(a < b && b < c)
    if InfoLevel>=10
        fprintf(' a=%-g \t b=%-g \t c=%-g \n',a,b,c)
        fprintf('parabolamin: values not ordered correctly ')
        fprintf(' a=%-g \t b=%-g \t c=%-g fa=%-g fb=%-g fc=%-g \n ',a,b,c,fa,fb,fc)
    end
    status=1; xmin=NaN;
    return
end

if isinf(fa) && isinf(fb) && isinf(fc) 
    status=1;
    xmin=a/100;
    warning('parabolamin:NaN',' On input fa, fb and fc are all NaN ')
    fprintf(' a=%-g \t b=%-g \t c=%-g fa=%-g fb=%-g fc=%-g \n ',a,b,c,fa,fb,fc)
    return
end

if isinf(fb) && isinf(fc)
    status=1;
    if ~isnan(b)
        xmin=b/100; % this suggests that the step size was too large so a starkly reduced one is return
    else
        xmin=1e-5;
    end
    warning('parabolamin:NaN',' On input fb and fc are NaN ')
    fprintf(' a=%-g \t b=%-g \t c=%-g fa=%-g fb=%-g fc=%-g \n ',a,b,c,fa,fb,fc)
    return
end


t1=(b-a)^2*(fb-fc);
t2=(b-c)^2*(fb-fa);

t3=(b-a)*(fb-fc);
t4=(b-c)*(fb-fa);

if t3 ~= t4
    xmin=b-0.5*(t1-t2)/(t3-t4);
else
    xmin=NaN; status=1;
end

%  rather than returning NaN I try to do something that might be reasonable in a one-sided line search
if isnan(xmin) || (fa+(b-a)*(fc-fa)/(c-a)<fb)  %  not convex
    
    if InfoLevel>=10000
        fprintf(' parabolic fit not successful, try to improve.\n  ')
    end
    if fc < fb && fb < fa      % decreasing values
        xmin=5*c;    % extrapolate
    elseif fc > fb && fb > fa  % increasing values
        xmin=a+b/10;
    elseif fc > fa && fb > fa  % fc=fb>fa
        xmin=a+b/100;
        fprintf(' xmin=a+b/100  \n')
    else                       %
        xmin=(a+b)/3;
    end
end


end

