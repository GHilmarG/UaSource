function Number=RoundNumber(Number,nDigits)


if nargin<2
    nDigits=2;
end

pm=sign(Number);
Number=abs(Number);

k=10^(floor(log10(Number)-nDigits+1));


Number=pm*k*round(Number/k);

end

 
 