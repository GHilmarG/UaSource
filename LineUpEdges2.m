function [xPolygon,yPolygon]=LineUpEdges2(CtrlVar,xa,xb,ya,yb,LineMax)


% [xPolygon,yPolygon]=LineUpEdges2(CtrlVar,xa,xb,ya,yb,LineMax)
%
% Lines up edges. 
% Takes edges (i.e. line segments) and lines them up to form continous lines, with NaN where there is a gap between edges
%
%
%  xa,  xb, ya, yb     : vectors defining the start and end x,y coordinates of the edges 
%                        for example: the n-th edge goes from [xa(n) ya(n)] to [xb(n) yb(n)]
%  LineMax             : maximum number of lined-up edges returned. Edges are returned in the order of the number of points in each edge.
%                        If, for example, LineMax=1, then only the single longest line is returned
%
%
% Note: As a part of the Mapping Toolbox there is a matlab routine `polymerge' that 
%       can also be used to do this, but this routine is
%       hoplessly slow and memory hungry for a large number of line segments.
%
%
%  
%
%  To plot:   plot(xPolygon,yPolygon)
%
%  To plot in different colours showing individual line segments
% figure ; 
% i=0;
% I=find(isnan(xPolygon)) ; 
% I=[1;I(:)];
% col=['b','r','c','g','k','m']; 
% for ii=1:numel(I)-1
%     i=i+1;
%     plot(xPolygon(I(ii):I(ii+1)),yPolygon(I(ii):I(ii+1)),col(i)) ; axis equal ; hold on ;
%     if i==numel(col) ; i=0 ; end
% end




Tolerance=100*eps;

if isempty(CtrlVar) 
    CtrlVar.InfoLevel=0;
end

if ~isfield(CtrlVar,'InfoLevel') ; CtrlVar.InfoLevel=0 ; end

if nargin<6
    LineMax=inf;
end

N=length(xa);
xPolygon=zeros(3*N,1)+inf ; yPolygon=zeros(3*N,1)+inf; % upper estimate, all edges are seperate lines segments

i=1; 
l=i; % l is the starting point of the current segment
xPolygon(i)=xa(1) ;  yPolygon(i)=ya(1) ;  i=i+1 ;
xPolygon(i)=xb(1) ;  yPolygon(i)=yb(1)  ; i=i+1 ;
xa(1)=NaN ; ya(1)=NaN ; xb(1)=NaN ; yb(1)=NaN ;


k=i-1;  % k is index of the current end of the merged line 
flipped=0;

while ~all(isnan(xa))
    
    
    % find input GL location closest to last included output GL pos
    % the current endpoint of the polygone is 
    % (xPolygon(k),yPolygon(k))
    sa=(xPolygon(k)-xa).^2+(yPolygon(k)-ya).^2;  % distance to all the `a' end points of remaining edges to the current end point of the merged line
    sb=(xPolygon(k)-xb).^2+(yPolygon(k)-yb).^2;  % distance to all the `b' end points of remaining edges to the current end point of the merged line
    [dista,ia]=min(sa);
    [distb,ib]=min(sb);
    
    %[xPolygon(:) yPolygon(:)]
    % if the distane is zero, then I am on the same GL
    % otherwise start a new GL
    if dista<= Tolerance
        %fprintf('dista==0\n')
        %xPolygon(i)=xa(ia) ;  yPolygon(i)=ya(ia) ;  i=i+1;
        xPolygon(i)=xb(ia) ;  yPolygon(i)=yb(ia) ;  i=i+1;
        k=i-1;
        xa(ia)=NaN ; ya(ia)=NaN ; xb(ia)=NaN ; yb(ia)=NaN ; % get rid of this edge as it now is a part of the merged line
    elseif distb<Tolerance
        %fprintf('distb==0\n')
        %xPolygon(i)=xb(ib) ;  yPolygon(i)=yb(ib) ;  i=i+1;
        xPolygon(i)=xa(ib) ;  yPolygon(i)=ya(ib) ;  i=i+1;
        k=i-1;
        xa(ib)=NaN ; ya(ib)=NaN ; xb(ib)=NaN ; yb(ib)=NaN ;
    elseif ~flipped
        
        flipped=1;
        % Now flip the line segment and transverse in the oposite direction
        xPolygon(l:k)=flipud(xPolygon(l:k));
        yPolygon(l:k)=flipud(yPolygon(l:k));
        
        
    
    else
        
        flipped=0;
        xPolygon(i)=NaN ;  yPolygon(i)=NaN ;  i=i+1; % NaN put between GLs
        
        l=i ; % l is the starting point of the current segment
        if dista<distb
            %fprintf('dista smaller\n')
            
            xPolygon(i)=xa(ia) ;  yPolygon(i)=ya(ia) ;  i=i+1;
            xPolygon(i)=xb(ia) ;  yPolygon(i)=yb(ia) ;  i=i+1;
            k=i-1;
            xa(ia)=NaN ; ya(ia)=NaN ; xb(ia)=NaN ; yb(ia)=NaN ;
        else
            %fprintf('distb smaller\n')
            
            xPolygon(i)=xb(ib) ;  yPolygon(i)=yb(ib) ;  i=i+1;
            xPolygon(i)=xa(ib) ;  yPolygon(i)=ya(ib) ;  i=i+1;
            k=i-1;
            xa(ib)=NaN ; ya(ib)=NaN ; xb(ib)=NaN ; yb(ib)=NaN ;
        end
    end
    
end

I=isinf(xPolygon);
xPolygon(I)=[] ; yPolygon(I)=[];

%% rearange GLs in order of total number of points withing each GL


xPolygon=[xPolygon;NaN] ; yPolygon=[yPolygon;NaN];
temp=sort(find(isnan(xPolygon))) ; temp=[0;temp;numel(xPolygon)];
[~,I]=sort(diff(temp),'descend');

NGL=min([numel(I),LineMax]);

if CtrlVar.InfoLevel>=10;
    fprintf('LineUpEdges: Found %-i grounding lines. Returning %-i.  \n',numel(I),NGL)
end
xx=[] ; yy=[]; 
for l=1:NGL
    
    n1=temp(I(l))+1 ; n2=temp(I(l)+1) ;
    xx=[xx;xPolygon(n1:n2)] ; yy=[yy;yPolygon(n1:n2)];
    
     
end

xx(end)=[] ;yy(end)=[];

xPolygon=xx; yPolygon=yy;






end
