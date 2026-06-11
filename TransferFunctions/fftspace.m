function freq=fftspace(n,dx) 

freq=zeros(n,1) ;
tpi=2*pi;
freq(1)=0.;
for i=2:n/2
    freq(i)=tpi*double(i-1)/double(n)/dx;
end
for i=n/2+1:n
    freq(i)=-tpi*double(n-i+1)/double(n)/dx;
end 

