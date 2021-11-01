%%
syms x y

% 3 node




N=3 ;

switch N
    
    case 3
        F{1}=x ;
        F{2}=1-x-y ;
        F{3}=y;
        
    case 6
        
        F{1}=(2*x-1)*x ;
        F{2}=4*x*(1-x-y) ;
        F{3}=(2*(1-x-y)-1)*(1-x-y) ;
        F{4}=4*y*(1-x-y);
        F{5}=(2*y-1)*y ;
        F{6}=4*x*y;
end

for I=1:N
    dF{I,1}=diff(F{I},x) ;
    dF{I,2}=diff(F{I},y) ;
end


for I=1:N
    ddF{I,1,1}=diff(dF{I,1},x) ;
    ddF{I,1,2}=diff(dF{I,1},y) ;
    ddF{I,2,1}=diff(dF{I,2},x) ;
    ddF{I,2,2}=diff(dF{I,2},y) ;
end