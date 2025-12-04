function [s]=FunctionSurface2d(Experiment,coordinates)
    
    
    
    x=coordinates(:,1) ; y=coordinates(:,2);
    Nnodes=length(x);
    
    switch Experiment
        case{'inverse1','forwardC','forwardb'}
            surfaceslope=0.01;
            s=-surfaceslope*(x-min(x(:)));
            
        otherwise
            disp(' case not recognised ')
    end
    
    
end