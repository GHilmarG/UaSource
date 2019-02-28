function sum=Looping(n,a,b)
%#codgen
%% codegen -config:mex Looping -args {1 , 1 , 1 }

sum=0;
sum2=0;
sum3=0;
for I=1:n
    
    a=a+1;
    %if sum>1e5
    %    sum=0;
    %end
    sum=sum+a*b;
    
    for j=1:n
        sum2=sum+b;
        
        for k=1:n
           sum3=sum3+sum2/2;
        end
    end
    
end


sum=sum+sum2+sum3;




return
end
