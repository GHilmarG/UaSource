function varargout=ExplicitEstimation(dt,dtRatio,Itime,varargin)

% varargout=ExplicitEstimation(dt,dtRatio,Itime,varargin)
% explicit estimate using second-order variable time step Adams-Bashforth
% http://lucan.me.jhu.edu/wiki/index.php/Second-order_variable_time_step_Adams-Bashforth_method
%


nInputs=nargin-3;
nOutputs = nargout;
varargout = cell(1,nOutputs);

if nInputs~=3*nOutputs
    error(' wrong number of inputs ')
end

if Itime==1 

    for I=1:nOutputs
        varargout{I}=varargin{1+3*(I-1)};
    end
    
elseif Itime==2 
    
    for I=1:nOutputs
        varargout{I}=varargin{1+3*(I-1)}+dt*varargin{2+3*(I-1)};
    end


elseif Itime >= 3  % second-order variable time step Adams-Bashforth
   
 
    
    
    for I=1:nOutputs
        varargout{I}=varargin{1+3*(I-1)}+dt*((1+0.5*dtRatio)*varargin{2+3*(I-1)}-0.5*dtRatio*varargin{3+3*(I-1)});
    end


end

end