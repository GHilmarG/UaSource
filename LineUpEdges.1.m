function [xEdge,yEdge,EdgeIndex]=LineUpEdges(xa,xb,ya,yb,LineMax)



% now just walk through the points always moving to the next closest point

N=length(xa);
xEdge=zeros(3*N,1)+Inf ; yEdge=zeros(3*N,1)+Inf; % upper estimate, all GL elements have seperate GLs

i=1;
xEdge(i)=xa(1) ;  yEdge(i)=yb(1) ;  i=i+1 ;
xEdge(i)=xa(1) ;  yEdge(i)=ya(1)  ; i=i+1 ;
xa(1)=[] ; ya(1)=[] ; xb(1)=[] ; yb(1)=[] ;


k=2; % last found output GL point

while ~isempty(xa)
    
    
    % find input GL location closest to last included output GL pos
    
    sa=(xEdge(k)-xa).^2+(yEdge(k)-ya).^2;
    sb=(xEdge(k)-xb).^2+(yEdge(k)-yb).^2;
    [dista,ia]=min(sa);
    [distb,ib]=min(sb);
    
    
    % if the distane is zero, then I am on the same GL
    % otherwise start a new GL
    if dista<= 100*eps
        %fprintf('dista==0\n')
        xEdge(i)=xa(ia) ;  yEdge(i)=ya(ia) ;  i=i+1;
        xEdge(i)=xb(ia) ;  yEdge(i)=yb(ia) ;  i=i+1;
        k=i-1;
        xa(ia)=[] ; ya(ia)=[] ; xb(ia)=[] ; yb(ia)=[] ;
    elseif distb<100*eps
        %fprintf('distb==0\n')
        xEdge(i)=xb(ib) ;  yEdge(i)=yb(ib) ;  i=i+1;
        xEdge(i)=xa(ib) ;  yEdge(i)=ya(ib) ;  i=i+1;
        k=i-1;
        xa(ib)=[] ; ya(ib)=[] ; xb(ib)=[] ; yb(ib)=[] ;
    else
        
        xEdge(i)=NaN ;  yEdge(i)=NaN ;  i=i+1; % NaN put between GLs
        
        if dista<distb
            %fprintf('dista smaller\n')
            
            xEdge(i)=xa(ia) ;  yEdge(i)=ya(ia) ;  i=i+1;
            xEdge(i)=xb(ia) ;  yEdge(i)=yb(ia) ;  i=i+1;
            k=i-1;
            xa(ia)=[] ; ya(ia)=[] ; xb(ia)=[] ; yb(ia)=[] ;
        else
            %fprintf('distb smaller\n')
            
            xEdge(i)=xb(ib) ;  yEdge(i)=yb(ib) ;  i=i+1;
            xEdge(i)=xa(ib) ;  yEdge(i)=ya(ib) ;  i=i+1;
            k=i-1;
            xa(ib)=[] ; ya(ib)=[] ; xb(ib)=[] ; yb(ib)=[] ;
        end
    end
    
end

I=isinf(xEdge);
xEdge(I)=[] ; yEdge(I)=[];

%% rearange GLs in order of total number of points withing each GL


xEdge=[xEdge;NaN] ; yEdge=[yEdge;NaN];
temp=sort(find(isnan(xEdge))) ; temp=[0;temp;numel(xEdge)];
[temp2,I]=sort(diff(temp),'descend');

NGL=min([numel(I),LineMax]);

fprintf('Found %-i grounding lines. Returning %-i.  \n',numel(I),NGL)

xx=[] ; yy=[]; EdgeIndex=1;
for l=1:NGL
    
    n1=temp(I(l))+1 ; n2=temp(I(l)+1) ;
    EdgeIndex=[EdgeIndex;numel(xx)+2+n2-n1];
    xx=[xx;xEdge(n1:n2)] ; yy=[yy;yEdge(n1:n2)];
    
end
EdgeIndex(end)=[]; xx(end)=[] ;yy(end)=[];

xEdge=xx; yEdge=yy;

%% There is a potential issue that I have not taken care of yet
% if  GL is not closed it is possible that it will be returned as two grounding lines..
%
%     if NGL>1
%         xStart=xGL(GLindex);
%         nGLend=[GLindex+2;numel(xGL)] ; nGLend(1)=[];
%         xEnd=xGL(nGLend);
%         [xStart xEnd]
%     end


end
