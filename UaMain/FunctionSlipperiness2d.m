function [C]=FunctionSlipperiness2d(Experiment,coordinates)
    
    
    
    x=coordinates(:,1) ; y=coordinates(:,2);
    Nnodes=length(x);
    
    switch Experiment
        case{'inverse1'}
            C=x*0+10;
        case{'forwardC'}
            C0=10; hmean=1000;
            ampl=0.1; sigma_x=5*hmean ; sigma_y=5*hmean ;
            DeltaC=ampl*exp(-((x/sigma_x).^2+(y/sigma_y).^2));
            DeltaC=DeltaC-mean(DeltaC);
            C=C0*(1+DeltaC);
        case{'forwardb'}
            C0=10; hmean=1000;
            ampl=0.0; sigma_x=5*hmean ; sigma_y=5*hmean ;
            DeltaC=ampl*exp(-((x/sigma_x).^2+(y/sigma_y).^2));
            DeltaC=DeltaC-mean(DeltaC);
            C=C0*(1+DeltaC);
            
        otherwise
            disp(' case not recognised ')
    end
    
    
end