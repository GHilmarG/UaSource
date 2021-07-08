
function [C,m,q,muk]=TestSlipperinessInputValues(CtrlVar,MUA,C,m,q,muk)

narginchk(4,6)
nargoutchk(2,4)



[nC,mC]=size(C);
[nm,mm]=size(m);


if nargin<5
    q=[];
end

if nargin<6
    muk=[];
end


if numel(C)==1
    %fprintf(' C given by user is a scalar. Assuming that C is same everywhere. \n')
    if  CtrlVar.CisElementBased
        C=C+zeros(MUA.Nele,1);
    else
        C=C+zeros(MUA.Nnodes,1);
    end
end


if numel(m)==1
    
    if  CtrlVar.CisElementBased
        m=m+zeros(MUA.Nele,1);
    else
        m=m+zeros(MUA.Nnodes,1);
    end
end



if isempty(q)
    
    pattern=["Budd","W-N0"];
    if contains(CtrlVar.SlidingLaw,pattern)
        fprintf("Input Error:  \t For sliding law: %s \n \t \t \t \t q must be defined in DefineSlipperiness.m \n",CtrlVar.SlidingLaw)
        fprintf("\t \t \t \t and in an inverse run in DefineInputsForInverseRun.m as well. \n")
        error("Incorrect inputs")
        
    end
    
else
    
    if numel(q)==1
        
        if  CtrlVar.CisElementBased
            q=q+zeros(MUA.Nele,1);
        else
            q=q+zeros(MUA.Nnodes,1);
        end
    end
    
end



if isempty(muk)
    
    pattern=["Tsai","Coulomb","Cornford","Umbi","minCW-N0","rpCW-N0","rCW-N0"]  ;
    if contains(CtrlVar.SlidingLaw,pattern)
        fprintf("Input Error:  \t For sliding law: %s \n \t \t \t \t muk must be defined in DefineSlipperiness.m \n",CtrlVar.SlidingLaw)
        fprintf("\t \t \t \t and in an inverse run in DefineInputsForInverseRun.m as well. \n")
        error("Incorrect inputs")
        
    end
    
else
    
    if numel(muk)==1
        
        if  CtrlVar.CisElementBased
            muk=muk+zeros(MUA.Nele,1);
        else
            muk=muk+zeros(MUA.Nnodes,1);
        end
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


[C,iU,iL]=kk_proj(C,CtrlVar.Cmax,CtrlVar.Cmin);

niU=numel(find(iU));
niL=numel(find(iL));

if CtrlVar.InfoLevel>=10
    if niU>0
        fprintf('TestCInputValues: On input %i C values are larger than largest allowed value (%g).  \n',niU);
        fprintf(' These values have been reset to the maximum allowed value of %g .\n',niU,CtrlVar.Cmax);
        
    end
    if niL>0
        fprintf('TestCInputValues: On input %i C values are smaller than smallest allowed value (%g). \n',niL,CtrlVar.Cmin);
        fprintf(' These values are have been reset to the minimum allowed value of %g.\n',CtrlVar.Cmin);
    end
end



end
