function [varargout]=InterpolateIntVariables(DTint,Iint,x,y,varargin)
    
    % interpolates integration points variables onto (x,y)
    
    Nin=4;  nIn=nargin-Nin;
    varargout=cell(nargout,1);
    
    
    for I=1:nIn
        fint=varargin{I};
        temp=fint(:) ; temp=temp(Iint);
        varargout{I}=Grid1toGrid2(DTint,temp,x,y);
    end
    
    
end