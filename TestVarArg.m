function [varargout]=TestVarArg(varargin)
    
    
    for I=1:nargin
        varargin{I}
    end
    varargout=cell(nargout,1);
    for I=1:nargout
        varargout{I}=I;
    end
    
end
