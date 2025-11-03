function [x,y,I] = Arrange2dPos(x,y)
    
    %
    % takes vectors of x and y values and tries to arrange them so that they form a line
    %
    
    if ~isempty(x) && ~isempty(y)
        
        xr=x ; yr=y;
        
        xc=mean(x) ; yc=mean(y);
        s=(x-xc).^2+(y-yc).^2;
        [~,iloc]=max(s);
        
        %[~,iloc]=min(y);
        
        N=numel(x);
        I=zeros(N,1);
        
        I(1)=iloc;
        xc=x(iloc) ; yc=y(iloc);
        x(iloc)=NaN; y(iloc)=NaN;
        
        % now just walk through the points always moving to the next closest point
        d=zeros(N-1,1);
        for k=1:N-1
            s=(x-xc).^2+(y-yc).^2;
            [dist,iloc]=min(s);
            d(k)=dist;
            I(k+1)=iloc;
            xc=x(iloc) ; yc=y(iloc);
            x(iloc)=NaN; y(iloc)=NaN; % don't reuse this point
        end
        

         x=xr(I) ; y=yr(I);
         
         % try to split GL up into individual GLs, and only return the longest one
         ind=find(d>10*std(d));
         
         
         if ~isempty(ind)
             x(ind+1)=NaN; y(ind+1)=NaN;
             
             temp=sort(find(isnan(x))) ; temp=[0;temp;numel(x)+1];
             
             [temp2,I]=max(diff(temp));
             
             n1=temp(I)+1 ;
             n2=temp(I+1)-1 ;
             x=x(n1:n2) ; y=y(n1:n2);
         end
    else
        %fprintf(' x and y empty on input to Arrange2dPos! \n')
        x=[] ; y=[] ; I=[] ;
    end
    
    %figure ; plot(x,y,'-o') ; axis equal
    %%
end
