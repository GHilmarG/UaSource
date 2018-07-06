function [coo,con] = NodalSweep(coo,con,alpha)
        
    
    temp=coo(:,1)*cos(alpha)+coo(:,2)*sin(alpha);
    [t,p]=sort(temp);
    
    p2=p*0; p2(p)=1:length(p2);
    coo(:,2)=coo(p,2);
    coo(:,1)=coo(p,1);
    con=p2(con);
            
    
end

