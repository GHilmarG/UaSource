function [coordinates,connectivity]=FEmeshSmoothing(coordinates,connectivity,maxit,tol)
    
%%
%
% Wrapper around the GHGsmoothmesh function, which in turn is a minor
% modification to a code by Darren Engwirda.
%
%
% Laplacian smoothin of the FE mesh.  This is an utility function.
%
%
%%

    [Nele,nod]=size(connectivity);
    
    [coordinates,connectivity]=ChangeElementType(coordinates,connectivity,3);
    
    % might introduce selected smoothing of elements here some day
    [coordinates] = GHGsmoothmesh(coordinates,connectivity,maxit,tol);
    
    [coordinates,connectivity]=ChangeElementType(coordinates,connectivity,nod);
    
end

