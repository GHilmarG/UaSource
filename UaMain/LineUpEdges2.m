function [xPolygon,yPolygon]=LineUpEdges2(CtrlVar,xa,xb,ya,yb,LineMax)

%% Line up line segments to form continous lines with NaN where gaps.
% [xPolygon,yPolygon]=LineUpEdges2(CtrlVar,xa,xb,ya,yb,LineMax)
%
% Lines up edges. Takes edges (i.e. line segments) and lines them up to form
% continous lines, with NaN where there is a gap between edges
%
%
%  xa,  xb, ya, yb     : vectors defining the start and end x,y coordinates of
%                        edges, for example: the n-th edge goes from [xa(n) ya(n)] to
%                        [xb(n) yb(n)]
%  LineMax             : maximum number of lined-up edges returned. Edges are
%                        returned in the order of the number of points in each edge.
%                        If, for example, LineMax=1, then only the single
%                        longest line is returned
%
%
% Line segments are considered to belong to seperate lines if the distance between start/end points
% is larger than CtrlVar.LineUpTolerance
%
% The default value is CtrlVar.LineUpTolerance=100*eps ;
%
%
% Note: As a part of the Mapping Toolbox there is a matlab routine `polymerge'
% that can also be used to do this, but this routine is hopelessly slow and
% memory hungry for a large number of line segments.
%
%  To plot:   plot(xPolygon,yPolygon)
%
%  To plot in different colours showing individual line segments
% figure ; i=0; I=find(isnan(xPolygon)) ; I=[1;I(:)];
% col=['b','r','c','g','k','m']; for ii=1:numel(I)-1
%     i=i+1; plot(xPolygon(I(ii):I(ii+1)),yPolygon(I(ii):I(ii+1)),col(i)) ; axis
%     equal ; hold on ; if i==numel(col) ; i=0 ; end
% end

%save TestSave ; error('fsda')

if numel(xa) ~= numel(xb)
    error('Ua:LineUpEdges:IncorrectNumberOfElements','number of elements in xa must equal number of elements in xb')
end

if numel(xa) ~= numel(ya)
    error('Ua:LineUpEdges:IncorrectNumberOfElements','number of elements in xa must equal number of elements in ya')
end


if numel(ya) ~= numel(yb)
    error('Ua:LineUpEdges:IncorrectNumberOfElements','number of elements in ya must equal number of elements in yb')
end


if isempty(xa)
    xPolygon=[] ; yPolygon=[];
    return
end




if isempty(CtrlVar)
    CtrlVar.InfoLevel=0;
    CtrlVar.LineUpTolerance=1000*eps ;
    
end

if ~isfield(CtrlVar,'InfoLevel') ; CtrlVar.InfoLevel=0 ; end
if ~isfield(CtrlVar,'LineUpTolerance') ; CtrlVar.LineUpTolerance=1000*eps ; end

% Tolerance controls how close endpoints of edges must be for them to be considered to be a part of the
% same polygon.
% Tolerance2 is a numerical tolerance used to determine if a x,y coordinate is identical to another
% x,y coordinate.

Tolerance=CtrlVar.LineUpTolerance;
Tolerance2=100*eps;

if nargin<6
    LineMax=inf;
end

% sort values in a circular manner around the mean centerpoint.
% This usually speeds things up as it makes it more likely that individual line segments are already aligned.
Theta=atan2(ya-mean(ya),xa-mean(xa));
[~,I]=sort(Theta);
xa=xa(I) ; ya=ya(I) ; xb=xb(I) ; yb=yb(I);

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
    [dista,ia]=min(sa);   % this ignores all NaN (of which there could be many)
    [distb,ib]=min(sb);
    
    %[xPolygon(:) yPolygon(:)]
    % if the distance is zero, then I am on the same GL
    % otherwise start a new GL
    
    %
    % Algorithm:
    % -Take edge points as defined by [xa ya ; xb yb] and add to Polygon
    % Find endpoints of edges closest to last added point in Polygon. If distance to either end points of
    % some edges is within a given tolerance, add those to Polygon and delete the edge.
    % If no such point is found, consider the possibility that there are some such edges close to the
    % other end of current Polygon and therefore flip the polygon once.
    % If no such edges are found, after having flipped once, start a new polygon, and define a new starting point.
    %
    if dista< Tolerance || distb<Tolerance
        
        if dista<= distb
            %fprintf('dista %f \n',dista)
            if dista>Tolerance2
                xPolygon(i)=xa(ia) ;  yPolygon(i)=ya(ia) ;  i=i+1;
            end
            xPolygon(i)=xb(ia) ;  yPolygon(i)=yb(ia) ;  i=i+1;
            k=i-1;
            xa(ia)=NaN ; ya(ia)=NaN ; xb(ia)=NaN ; yb(ia)=NaN ; % get rid of this edge as it now is a part of the merged line
            
            while ia<=(numel(xa)-1)
                %aCase=aCase+1;
                Test=abs(xa(ia+1)-xPolygon(k))+abs(ya(ia+1)-yPolygon(k));
                if Test<Tolerance2
                    ia=ia+1;
                    xPolygon(i)=xb(ia) ;  yPolygon(i)=yb(ia) ;  i=i+1;
                    k=i-1;
                    xa(ia)=NaN ; ya(ia)=NaN ; xb(ia)=NaN ; yb(ia)=NaN ;
                else
                    break
                end
            end
        else
            %fprintf('distb==0\n')
            %xPolygon(i)=xb(ib) ;  yPolygon(i)=yb(ib) ;  i=i+1;
            
            if dista>Tolerance2
                xPolygon(i)=xb(ib) ;  yPolygon(i)=yb(ib) ;  i=i+1;
            end
            
            xPolygon(i)=xa(ib) ;  yPolygon(i)=ya(ib) ;  i=i+1;
            k=i-1;
            xa(ib)=NaN ; ya(ib)=NaN ; xb(ib)=NaN ; yb(ib)=NaN ;
            
            while ib<=(numel(xa)-1)
                %bCase=bCase+1;
                Test=abs(xa(ib+1)-xPolygon(k))+abs(ya(ib+1)-yPolygon(k));
                if Test<Tolerance2
                    ib=ib+1;
                    xPolygon(i)=xa(ib) ;  yPolygon(i)=ya(ib) ;  i=i+1;
                    k=i-1;
                    xa(ib)=NaN ; ya(ib)=NaN ; xb(ib)=NaN ; yb(ib)=NaN ;
                else
                    break
                end
            end
        end
        
    elseif ~flipped   % flipp once
        
        flipped=1;
        % Now flip the line segment and transverse in the oposite direction
        xPolygon(l:k)=flipud(xPolygon(l:k));
        yPolygon(l:k)=flipud(yPolygon(l:k));
        
        
        
    else  % dist larger than tolerance, must start a new line-segment.
        % Put NaN to mark the division between line segments, and find a starting point
        % for the next line segment.
        
        
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
    
    %     figure(10) ; hold on  ; plot(xPolygon,yPolygon,'o-')
    %     prompt = 'Do you want more? Y/N [Y]: ';
    %     [xa(:) ya(:) xb(:) yb(:)]
    %     str = input(prompt,'s');
    %     if isempty(str)
    %         str = 'Y';
    %     end
    
end

I=isinf(xPolygon);
xPolygon(I)=[] ; yPolygon(I)=[];

%% rearange GLs in order of total number of points withing each GL


xPolygon=[xPolygon;NaN] ; yPolygon=[yPolygon;NaN];
temp=sort(find(isnan(xPolygon))) ; temp=[0;temp;numel(xPolygon)];
[~,I]=sort(diff(temp),'descend');

NGL=min([numel(I),LineMax]);

if CtrlVar.InfoLevel>=10
    fprintf('LineUpEdges2: Found %-i seperate lines. Returning %-i.  \n',numel(I),NGL)
end


M=10000;
xx=zeros(M,1) ; yy=zeros(M,1) ;
N=0;
for l=1:NGL
    
    n1=temp(I(l))+1 ;
    n2=temp(I(l)+1) ;
    n=n2-n1+1;
    
    
    if length(xx)<(N+n)
        M=2*M;
        xx=[xx;zeros(M,1)];
        yy=[yy;zeros(M,1)];
    end
    
    xx(N+1:N+n)=xPolygon(n1:n2) ; yy(N+1:N+n)=yPolygon(n1:n2) ;
    N=N+n;
    
    
end


xx(N:end)=[] ;  yy(N:end)=[] ;



%% Added in 2022 January, extra checks to get rid of short line segments and segments with just one point

% Now get rid of very short line segments:
ds = sqrt(sum(diff([xx yy],[],1).^2,2)) ;
I=ds<CtrlVar.LineUpTolerance ;
xx(I)=[] ; yy(I)=[] ;

% And now get rid of any single -point "line" segments:
I=isnan(xx);
K=find(I);                % labels of nan
if any(diff(K)==2)  % OK so there are some segments with just one point

    d=xx+nan ;                % generous pre-allocation
    nn=1;
    for k=1:numel(K)-1      % loop over nan
        iD=K(k+1)-K(k) ;    % if two adjacent nan are diff by 2, then the line between those has only one point

        % d is an array marking the points in xx and yy for deletion
        if iD==2
            d(nn)=K(k)+1;    % get rid of the data point
            d(nn+1)=K(k+1) ; % and one of the nan bracketing that data point
            nn=nn+2;
        end
    end

    % Is there a single data point at the end?  The above approach would not find that case
    if isnan(xx(end-1)) && ~isnan(xx(end))
        d(nn)=numel(xx)-1;
        nn=nn+1;
        d(nn)=numel(xx);
    end
    d(isnan(d))=[] ;

    xx(d)=[] ; yy(d)=[];  % and now delete from xx and yy



end


if isnan(xx(end))
    xx(end)=[];
    yy(end)=[];
end

%%

xPolygon=xx; yPolygon=yy;


end
