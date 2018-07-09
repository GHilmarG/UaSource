
function [C,m]=TestSlipperinessInputValues(CtrlVar,MUA,C,m)

[nC,mC]=size(C);
[nm,mm]=size(m);



if numel(C)==1
    %fprintf(' C given by user is a scalar. Assuming that C is same everywhere. \n')
    if  CtrlVar.CisElementBased
        C=C+zeros(MUA.Nele,1);
    else
        C=C+zeros(MUA.Nnodes,1);
    end
end


if numel(m)==1
    %fprintf(' C given by user is a scalar. Assuming that C is same everywhere. \n')
    if  CtrlVar.CisElementBased
        m=m+zeros(MUA.Nele,1);
    else
        m=m+zeros(MUA.Nnodes,1);
    end
end



if CtrlVar.AutomaticallyMapAGlenBetweenNodesAndEleIfEnteredIncorrectly
    
    if CtrlVar.CisElementBased
        
        if nC==MUA.Nnodes
            fprintf('\n Note:  C is element based, but entered on input as a nodal variable.\n')
            fprintf('        C will be mapped from nodes to elements by averaging over neighbouring nodes. \n')
            C=Nodes2EleMean(MUA.connectivity,C);
        end
        
        if nm==MUA.Nnodes
            fprintf('\n Note:  m is element based, but entered on input as a nodal variable.\n')
            fprintf('        m will be mapped from nodes to elements by averaging over neighbouring nodes. \n')
            m=Nodes2EleMean(MUA.connectivity,m);
        end
        
    else
        
        if nC==MUA.Nele || nm==MUA.Nele
            
            M=Ele2Nodes(MUA.connectivity,MUA.Nnodes);
            
            if nC==MUA.Nele
                
                fprintf('\n Note:  C is nodal based, but entered on input as an element variable.\n')
                fprintf('        C will be mapped from elements to nodes by averaging over neighbouring elements. \n')
                
                C=M*C;
            end
            
            if nm==MUA.Nele
                
                fprintf('\n Note:  m is nodal based, but entered on input as an element variable.\n')
                fprintf('        m will be mapped from elements to nodes by averaging over neighbouring elements. \n')
                
                m=M*m;
            end
        end
    end
end




if CtrlVar.CisElementBased  && ~(length(MUA.connectivity)==length(C))
    save TestSave ;
    error(' C is element-based but on input does not have same number of elements as there are elements in mesh. All variables saved in TestSave.mat ')
elseif ~CtrlVar.CisElementBased && ~(length(MUA.coordinates) == length(C))
    save TestSave ;
    error(' C is node-based but input does not have same number of elements as there are nodes in mesh. All variables saved in TestSave.mat ')
    
end


if CtrlVar.CisElementBased  && ~(length(MUA.connectivity)==length(m))
    save TestSave ;
    error('The variable m is element-based but on input does not have same number of elements as there are elements in mesh. All variables saved in TestSave.mat ')
elseif ~CtrlVar.CisElementBased && ~(length(MUA.coordinates) == length(m))
    save TestSave ;
    error('The variable m is node-based but input does not have same number of elements as there are nodes in mesh. All variables saved in TestSave.mat ')
end





end
