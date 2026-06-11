function rfrac=ReminderFraction(x,y)
    
    % ReminderFraction(x,y)
    % returns the ratio abs(x-y*round(x/y))/y
    % good for estimating if x is close to being a positive or negative integer multiple of y
    % 
    % if y->0 then all values become arbitrarily close to being an interger mupliple of y
    % if on input y<0 then rfrac=NaN

    
    rfrac=abs(x-y.*round(x./y))./y ;
    rfrac(y<=0)=NaN;

end
