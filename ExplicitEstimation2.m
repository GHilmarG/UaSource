function varargout=ExplicitEstimation2(dt,dtRatio,Itime,extrapolation,varargin)

% varargout=ExplicitEstimation(dt,dtRatio,Itime,varargin)
% explicit estimate using second-order variable time step Adams-Bashforth
% http://lucan.me.jhu.edu/wiki/index.php/Second-order_variable_time_step_Adams-Bashforth_method
%


if isempty(extrapolation)
    extrapolation='auto';
end


if isequal(extrapolation,'auto')
    if Itime==1
        extrapolation='constant';
    elseif Itime==2
        extrapolation='linear';
    elseif Itime >= 3  % second-order variable time step Adams-Bashforth
        extrapolation='quadratic';
    end
else
    if Itime==2 && isequal(extrapolation,'quadratic')
        extrapolation='linear';
    end
end

if Itime==1
    extrapolation='constant';
end

nInputs=nargin-4;
nOutputs = nargout;
varargout = cell(1,nOutputs);

if nInputs~=3*nOutputs
    error(' wrong number of inputs ')
end

switch extrapolation
    
    case 'constant'
        
        
        for I=1:nOutputs
            varargout{I}=varargin{1+3*(I-1)};
        end
        
    case 'linear'
        
        for I=1:nOutputs
            varargout{I}=varargin{1+3*(I-1)}+dt*varargin{2+3*(I-1)};
        end
        
    case 'quadratic'
        
        for I=1:nOutputs
            varargout{I}=varargin{1+3*(I-1)}+dt*((1+0.5*dtRatio)*varargin{2+3*(I-1)}-0.5*dtRatio*varargin{3+3*(I-1)});
        end
        
        
end

end