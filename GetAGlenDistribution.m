function [AGlen,n]=GetAGlenDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

[AGlen,n]=DefineAGlenDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);


if numel(AGlen)==1
    
    fprintf(' AGlen given by user is a scalar. Assuming that AGlen is same everywhere. \n')
    if  CtrlVar.AGlenisElementBased
        AGlen=AGlen+zeros(MUA.Nele,1);
    else
        AGlen=AGlen+zeros(MUA.Nnodes,1);
    end
    
end

if CtrlVar.AGlenisElementBased  && ~(length(MUA.connectivity)==length(AGlen))
    save TestSave ;
    error(' AGlen is element-based but does not have same number of elements as there are elements in mesh. All variables saved in TestSave.mat ')
elseif ~CtrlVar.AGlenisElementBased && ~(length(MUA.coordinates) == length(AGlen))
    save TestSave ;
    error(' AGlen is node-based but does not have same number of elements as there are nodes in mesh. All variables saved in TestSAve.mat ')
end

[AGlen,iU,iL]=kk_proj(AGlen,CtrlVar.AGlenmax,CtrlVar.AGlenmin);


end