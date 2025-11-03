


%[xPolygon,yPolygon]=LineUpEdges2(CtrlVar,xa,xb,ya,yb,LineMax)

%% not working! use LineUpEdges2  
Points=[0 0 ; 1 0 ; 2 0 ; 2 1; 20 0 ; 21 0];

xa=Points(1:end-1,1) ;  ya=Points(1:end-1,2) ; 
xb=Points(2:end,1) ;  yb=Points(2:end,2) ; 


points=[ [xa;xb] [ya;yb] ]; 

points=unique(points,'rows');

%Mdl=KDTreeSearcher(points);


[Idx, D] = knnsearch(Points,Points,'k',3);

Tolerance=10;

N=size(Points,1);
I=zeros(N,1)+NaN;
Used=false(N,1);

m=1;
I(m)=1; 
Used(m)=true;
m=m+1;

for k=1:N
        
    if ~Used(Idx(k,2)) && D(k,2) < Tolerance
       I(m)=Idx(k,2); 
       Used(I(m))=true;
       m=m+1;
    elseif ~Used(Idx(k,3)) && D(k,3) < Tolerance
        I(m)=Idx(k,3); 
        Used(k)=true; 
        m=m+1;
    elseif D(k,2)< Tolerance && D(k,3)< Tolerance
        % end of line segment
        I(m)=NaN ; 
        m=m+1;
        I(m)=k;
    else
        I(m)=k; 
        Used(k)=true; 
        m=m+1;
    end
        
    
end

Idx
I