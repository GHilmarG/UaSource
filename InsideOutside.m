function I=InsideOutside(xy,boundary)

%%
% Determines if xy is inside or outside of boundary
%
% Just a simple wrapper around inpoly2 to take care of the possibility of
% several boundaries defined by boundary, seperated by NaNs
%
%%

K=find(isnan(boundary(:,1)));


if numel(K)==0
    
    I=inpoly2(xy,boundary);
    
else
    
    I=false(size(xy(:,1)));
    i1=1;
    
    for j=1:numel(K)
        
        i2=K(j)-1;
        J=inpoly2(xy,boundary(i1:i2,:));
        i1=i2+2;
        
        I=I | J;
    end
    
    J=inpoly2(xy,boundary(i1:end,:));
    I=I | J;
end


end