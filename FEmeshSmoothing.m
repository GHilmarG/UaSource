function [coordinates,connectivity]=FEmeshSmoothing(coordinates,connectivity,maxit,tol)
    
    [Nele,nod]=size(connectivity);
    
    [coordinates,connectivity]=ChangeElementType(coordinates,connectivity,3);
    
    % might introduce selected smoothing of elements here some day
    [coordinates] = GHGsmoothmesh(coordinates,connectivity,maxit,tol);
    
    [coordinates,connectivity]=ChangeElementType(coordinates,connectivity,nod);
    
end

