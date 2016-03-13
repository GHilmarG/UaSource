function [C,m]=GetSlipperyDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

[C,m]=DefineSlipperyDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);



if numel(C)==1
    
    fprintf(' C given by user is a scalar. Assuming that C is same everywhere. \n')
    if  CtrlVar.CisElementBased
        C=C+zeros(MUA.Nele,1);
    else
        C=C+zeros(MUA.Nnodes,1);
    end
    
end

[C,iU,iL]=kk_proj(C,CtrlVar.Cmax,CtrlVar.Cmin);


end




