function [T,trep,WS,stats]=tsearchGHG(x,y,TRI,xi,yi,trep,WS)
    
    % [T,trep,WS,stats]=tsearch(x,y,TRI,xi,yi,trep,WS)
    % replacement for matlab tsearch
    % both trep and WS are optional and can also be defined as empty
    % if the triangulation does not change then use trep and WS from a previous call in all subsequent calls (considerable speed up)
    persistent PointsSize TriSize Ptrep PWS
    
    fprintf(' tsearch replacement starts ') ;  tic;
    
    x=x(:) ; y=y(:);
    ELEMS = int32(TRI)';
    markers = [xi(:) yi(:)]';
    points=[x(:)' ; y(:)'];
    
    if nargin<6  || isempty(trep) || isempty(WS)
        if ~isempty(PWS) & (size(points)==PointsSize) & (size(TRI)==TriSize)
            trep=Ptrep ;
            WS=PWS;
        else
            
            trep = TriRep(TRI, x,y);
            WS = [];
            WS.NEIGHBORS = trep.neighbors()';
            WS.NEIGHBORS(isnan(WS.NEIGHBORS)) = 0;
            WS.NEIGHBORS = int32(WS.NEIGHBORS);
            WS.xmin = min(x);
            WS.xmax = max(x);
            WS.ymin = min(y);
            WS.ymax = max(y);
            
            Ptrep=trep ;
            PointsSize=size(points);
            TriSize=size(TRI);
        end
    else
        Ptrep=trep ;
        PointsSize=size(points);
        TriSize=size(TRI);
    end
    setenv('OMP_NUM_THREADS', '1');
    %setenv('OMP_NUM_THREADS', '2');
    
    [T , WS , stats] = tsearch2(points, ELEMS, markers, WS);
    
    PWS=WS;  % have to do this hier because WS constaines qtree on output
        
    % final modifications to make output identical to that of tsearch
    T=T(:);
    T=double(T);
    T(T==0)=NaN;
    tcpu=toc;
    fprintf(' and returns %-g sec later \n ',tcpu)
end

