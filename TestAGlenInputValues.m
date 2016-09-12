function [AGlen,n]=TestAGlenInputValues(CtrlVar,MUA,AGlen,n)

%%
% [AGlen,n]=TestAGlenInputValues(CtrlVar,AGlen,n)
%
% Performs simple test to see if user input for A and n are OK
%

[nA,mA]=size(AGlen);
[nn,mn]=size(n);


if numel(AGlen)==1 || numel(n)==1
    
    %fprintf(' AGlen given by user is a scalar. Assuming that AGlen is same everywhere. \n')
    if  CtrlVar.AGlenisElementBased
        AGlen=AGlen+zeros(MUA.Nele,1);
        n=n+zeros(MUA.Nele,1);
    else
        AGlen=AGlen+zeros(MUA.Nnodes,1);
        n=n+zeros(MUA.Nnodes,1);
    end
    
end

if CtrlVar.AutomaticallyMapAGlenBetweenNodesAndEleIfEnteredIncorrectly
    
    if CtrlVar.AGlenisElementBased
        
        if nA==MUA.Nnodes
            fprintf('\n Note:  AGlen is element based, but entered on input as a nodal variable.\n')
            fprintf('        AGlen will be mapped from nodes to elements by averaging over neighbouring nodes. \n')
            AGlen=Nodes2EleMean(MUA.connectivity,AGlen);
        end
        
        if nn==MUA.Nnodes
            fprintf('\n Note:  n is element based, but entered on input as a nodal variable.\n')
            fprintf('        n will be mapped from nodes to elements by averaging over neighbouring nodes. \n')
            n=Nodes2EleMean(MUA.connectivity,n);
        end
        
    else
        
        if nA==MUA.Nele || nn==MUA.Nele
            
            M=Ele2Nodes(MUA.connectivity,MUA.Nnodes);
            
            if nA==MUA.Nele
                
                fprintf('\n Note:  AGlen is nodal based, but entered on input as an element variable.\n')
                fprintf('        AGlen will be mapped from elements to nodes by averaging over neighbouring elements. \n')
                
                AGlen=M*AGlen;
            end
            
            if nn==MUA.Nele
                
                fprintf('\n Note:  n is nodal based, but entered on input as an element variable.\n')
                fprintf('        n will be mapped from elements to nodes by averaging over neighbouring elements. \n')
                
                n=M*n;
            end
        end
    end
end


if CtrlVar.AGlenisElementBased  && ~(length(MUA.connectivity)==length(AGlen))
    save TestSave ;
    error(' AGlen is element-based but does not have same number of elements as there are elements in mesh. All variables saved in TestSave.mat ')
elseif ~CtrlVar.AGlenisElementBased && ~(length(MUA.coordinates) == length(AGlen))
    save TestSave ;
    error(' AGlen is node-based but does not have same number of elements as there are nodes in mesh. All variables saved in TestSave.mat ')
end

[AGlen,iU,iL]=kk_proj(AGlen,CtrlVar.AGlenmax,CtrlVar.AGlenmin);


niU=numel(find(iU));
niL=numel(find(iL));

if niU>0
    fprintf(' On input %i AGlen values are larger than largest allowed value (%g%). These values are reset.\n',niU,CtrlVar.AGlenmax);
end


if niL>0
    fprintf(' On input %i AGlen values are smaller than smallest allowed value (%g%). These values are reset.\n',niL,CtrlVar.AGlenmin);
end











end