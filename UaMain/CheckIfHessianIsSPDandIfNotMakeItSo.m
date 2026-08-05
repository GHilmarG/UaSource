


function [HlE,lEnd]=CheckIfHessianIsSPDandIfNotMakeItSo(H,MUA,lStart)


%%
%
% Can I use algorithm 4.3 to find lambda?
%
% https://digital.library.unt.edu/ark:/67531/metadc283525/m2/1/high_res_d/metadc283525.pdf
%

l=lStart;
lmin=lStart/1e6;


% First check the special case that lStart is zero, and H is pos def.
if l==0
    [~, flag] = chol(H);
    if flag==0

        fprintf("H is pos def, i.e.  for l=0 \n")
        HlE=H ;
        lEnd=0;
        return
    end

end


nNodes=size(MUA.M,1);
nH=size(H,1);

if nH==nNodes
    E=MUA.M ;
elseif nH==2*nNodes
    E=blkdiag(MUA.M,MUA.M) ;
else
    error("wrong dimentions")
end

% e=eigs(H) ;
% FindOrCreateFigure("Hessian eigenvalues") ; plot(sort(e),".r")
%

% check if H +l E is positive definite for the input value of l







if l==0
    l=1e-15;
end

HlE = H + l * E;
[~, flag] = chol(HlE);
 factor=100;

if flag==0 % H + l E is positive definite, but can I reduce l?

    while ~flag  && l > lmin    % decreasing case

        l=l/factor;
        HlE = H + l * E;
        [~, flag] = chol(HlE);
        fprintf(sprintf("Decreasing l until no longer positive definite. l=%g \n",l));

    end

    lEnd=l*factor;  % this was the value before it failed.
    fprintf("Hessian is positive definite for l=%g \n",l);

else


    while true  % H = L E was not pos, so I need to increase l : increasing case

        l=l*10;
        HlE = H + l * E;
        [~, flag] = chol(HlE);

        if flag == 0
            fprintf("Hessian is positive definite for l=%g \n",l);
            break;
        else
            fprintf(sprintf("Modifying Hessian to make it positive definite. l=%g \n",l));
        end
    end

    lEnd=l;

end


return

end

