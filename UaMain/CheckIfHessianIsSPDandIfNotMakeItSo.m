


function [HlE,lEnd]=CheckIfHessianIsSPDandIfNotMakeItSo(CtrlVar,H,E,lStart)

narginchk(4,4)

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
        if CtrlVar.InfoLevelInverse>=10
            fprintf("H is pos def, i.e.  for l=0 \n")
        end
        HlE=H ;
        lEnd=0;
        return
    end

end



if l==0
    l=1e-15;
end

HlE = H + l * E;
[~, flag] = chol(HlE);
DownFactor=10;
UpFactor=100;
iDecreaseMax=1;  % only decrease l by the fraction factor^iDecreaseMax each time
iDecrease=0;
if flag==0 % H + l E is positive definite, but can I reduce l?

    while ~flag  && l > lmin   && iDecrease<iDecreaseMax  % decreasing case

        l=l/DownFactor;
        HlE = H + l * E;
        [~, flag] = chol(HlE);
        iDecrease=iDecrease+1;
        if CtrlVar.InfoLevelInverse>=10

            fprintf(sprintf("%i: Decreasing l until no longer positive definite. l=%g \n",iDecrease,l));
        end

    end

    if flag~=0
        lEnd=l*DownFactor;  % this was the value before it failed.
        if CtrlVar.InfoLevelInverse>=10
            fprintf("Hessian is positive definite for l=%g \n",l);
        end
    else
        lEnd=l;
    end

else


    while true  % H = L E was not pos, so I need to increase l : increasing case

        l=l*UpFactor;
        HlE = H + l * E;
        [~, flag] = chol(HlE);

        if flag == 0
            if CtrlVar.InfoLevelInverse>=10
                fprintf("Hessian is positive definite for l=%g \n",l);
            end
            break;
        else
            if CtrlVar.InfoLevelInverse >=10
                fprintf(sprintf("Modifying Hessian to make it positive definite. l=%g \n",l));
            end
        end
    end

    lEnd=l;

end


return

end

