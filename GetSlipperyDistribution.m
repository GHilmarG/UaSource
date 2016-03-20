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


if CtrlVar.CisElementBased  && ~(length(MUA.connectivity)==length(C))
    save TestSave ;
    error(' C is element-based but on input does not have same number of elements as there are elements in mesh. All variables saved in TestSave.mat ')
elseif ~CtrlVar.CisElementBased && ~(length(MUA.coordinates) == length(C))
    save TestSave ;
    error(' C is node-based but input does not have same number of elements as there are nodes in mesh. All variables saved in TestSave.mat ')
    
end


[C,iU,iL]=kk_proj(C,CtrlVar.Cmax,CtrlVar.Cmin);


end




