function [varargout]=InterpolateNodalVariables(DTxy,x,y,varargin)
    
    % interpolates nodal variables onto (x,y)
    N=size(varargin,2);

    
    varargout=cell(1,N);
    
    
    for I=1:N
        varargout{I}=Grid1toGrid2(DTxy,varargin{I},x,y);
    end
    
    
end