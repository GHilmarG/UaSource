


function [H,Hcond]=ImproveMatrixCondition(H,E,l,lMin)

narginchk(4,4)


Hcond=condest(H);

if Hcond < 1e12

    if lMin==0

        return
    else

        l=lMin;  % this is now the l used, and it is used even if the condition was OK

    end


end

l=max(l,lMin);

nNodes=size(E,1);
nH=size(H,1);

if nH==nNodes

elseif nH==2*nNodes
    E=blkdiag(E,E) ;
else
    error("wrong dimentions")
end


scale=mean(abs(diag(H)))/ mean(abs(diag(E)));

H=H+l*scale*E;

Hcond2=condest(H);

fprintf("conditions improved from %g to %g or by the ratio of %g \n ",Hcond,Hcond2,Hcond/Hcond2)

Hcond=Hcond2;

return

end

