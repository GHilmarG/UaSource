function ticks=logticks(x,nPowRange,nTicks)

% just a rough way of generating approximatly sensible tickmarks for a log scale
% on input x is a set of numbers
%
%
%
%save TestSave

if min(x)<0
    ticks=[];
    return
end

if nargin<2 
    nPowRange=10;
end

if nargin<3
    nTicks=12;
end

if all(x==0)
    MinExp=-1 ; MaxExp=1;
else
    
    MinExp=floor(log10(min(x)));
    MaxExp=ceil(log10(max(x)));
end


if MaxExp-MinExp > nPowRange  % do not cover more than nPowRange orders of magnitudes
    MinExp=MaxExp-nPowRange;
end

Ticks=logspace(MinExp,MaxExp,MaxExp-MinExp+1);  % tickmarks at every power of 10
ticks=Ticks;



tv=[5 2 3 1.25 6 4 8 7 1.5 9 1.75]; iCount=0 ;
while numel(ticks)< nTicks && iCount<numel(tv)
    iCount=iCount+1;
    temp=Ticks*tv(iCount) ; temp=temp(1:end-1) ; ticks=sort([ticks temp]);
    
    if numel(ticks)>2
        if ticks(end-1)>=max(x) ; ticks=ticks(1:end-1) ; end
    end
end


end