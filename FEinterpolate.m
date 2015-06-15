function [h,s,b,B,a,u,v,etaInt,C]=FEinterpolate(Experiment,connectivity,coordinates,x,y,xint,yint,s,b,B,a,u,v,etaInt,C,nip);
    
    
    % take s, b, B and a defined on x y and interpolate on the new grid as defined by connectivity and coordinates
    
    xnew=coordinates(:,1) ; ynew=coordinates(:,2);
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    
    DT = DelaunayTri(x,y);  TRI=DT.Triangulation;
    
    
    F = TriScatteredInterp(DT,s,'natural'); s = F(xnew,ynew);
    F = TriScatteredInterp(DT,a,'natural'); a = F(xnew,ynew);
    F = TriScatteredInterp(DT,C,'natural'); C = F(xnew,ynew);
    F = TriScatteredInterp(DT,u,'natural'); u = F(xnew,ynew);
    F = TriScatteredInterp(DT,v,'natural'); v = F(xnew,ynew);
    
    B=FunctionBedrock2d(Experiment,coordinates);
    
    if strcmp(Experiment,'ex1a')
       F = TriScatteredInterp(DT,b,'natural'); b = F(xnew,ynew);
    else
       b=B; 
    end
       
    h=s-b;
    
    
    [xintNew,yintNew] = CalculateIntegrationPointCoordinatesVector(connectivity,coordinates,nip);
    
    meanEtaInt=mean(etaInt(:));
    F = TriScatteredInterp(xint(:),yint(:),etaInt(:)); etaInt = F(xintNew,yintNew);
    
    etaInt=reshape(etaInt,Nele,nip);
    
    % sometimes due to rounding errors the integration points are outside of the convex hull
    % This seems to happen if nod=3 and nip=3, where the integraton points are along the edges of the
    % triangles
    if any(isnan(b)) ; error(' FEinterpolate: interpolation gave rise to nan in b');  end
    if any(isnan(s)) ; error(' FEinterpolate: interpolation gave rise to nan in s');  end
    if any(isnan(a)) ; error(' FEinterpolate: interpolation gave rise to nan in a');  end
    if any(isnan(C)) ; error(' FEinterpolate: interpolation gave rise to nan in C');  end
    if any(isnan(u)) ; error(' FEinterpolate: interpolation gave rise to nan in u');  end
    if any(isnan(v)) ; error(' FEinterpolate: interpolation gave rise to nan in v');  end
    
    if any(isnan(etaInt))
        disp(' FEinterpolate: interpolation gave rise to nan in eta Int')
        etaInt(isnan(etaInt))=meanEtaInt;
    end
    
    
    
end


    
    
    
    
    