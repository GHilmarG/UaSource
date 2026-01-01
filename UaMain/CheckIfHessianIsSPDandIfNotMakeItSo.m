


function [H,lEnd]=CheckIfHessianIsSPDandIfNotMakeItSo(H,MUA,lStart)

lEnd=0;
%return


lStart=1e-8;
l=lStart;


nNodes=size(MUA.M,1);
nH=size(H,1);

if nH==nNodes
    EYE=MUA.M ;
elseif nH==2*nNodes
    EYE=blkdiag(MUA.M,MUA.M) ;
else
    error("wrong dimentions")
end


while true
    [~, pp] = chol(H);
    if pp == 0
        fprintf('Hessian is positive definite.\n');
        l=l/10;
        break;
    else
        fprintf('Modifying Hessian to make it positive definite.\n');
        H = H + l * EYE;
        l = l * 10;
    end
end

lEnd=l;

return

end

