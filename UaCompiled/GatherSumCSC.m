function vs=GatherSumCSC(Xval,invPtr,invIdx)
% vs(k) = sum of Xval(invIdx(invPtr(k)+1 : invPtr(k+1))), k=1..numel(invPtr)-1, in the order given by invIdx (see uvhCSCPattern.m).
% Compiled with MATLAB Coder (GatherSumCSC_mex), parallel over k with OpenMP: each k is written by one thread only, so there are no
% write conflicts, and the result does not depend on the number of threads.
nk=numel(invPtr)-1;
vs=zeros(nk,1);
parfor k=1:nk
    s=0;
    for p=invPtr(k)+1:invPtr(k+1)
        s=s+Xval(invIdx(p));
    end
    vs(k)=s;
end
end
