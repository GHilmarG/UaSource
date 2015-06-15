
function [s,b,S,B,alpha]=DefineGeometry(Experiment,CtrlVar,MUA,time,FieldsToBeDefined)
    % FieldsToBeDefined='sbSB' ; 
    
    x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
    
    alpha=0.01; hmean=1000; 
    ampl_b=0.5*hmean; sigma_bx=5000 ; sigma_by=5000;
    Deltab=ampl_b*exp(-((x/sigma_bx).^2+(y/sigma_by).^2));
    Deltab=Deltab-mean(Deltab);
    
    B=zeros(MUA.Nnodes,1) + Deltab;
    S=B*0-1e10;
    b=B;
    s=B*0+hmean;
    
    
end
