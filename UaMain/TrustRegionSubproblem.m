

function l=TrustRegionSubproblem(H,E,g,l,D)


if isempty(E)

    nNodes=size(MUA.M,1);
    nH=size(H,1);

    if nH==nNodes
        E=MUA.M ;
    elseif nH==2*nNodes
        E=blkdiag(MUA.M,MUA.M) ;
    else
        error("wrong dimentions")
    end

end



for iLoop=1:5

    HlE = H + l * E;

    [R, flag] = chol(HlE);  % change this later to use permutation matrix P


    if flag~=0
        fprintf("H+ l E not pos definite.\n")
    end

    p=R\(R'\(-g)) ;
    q=R'\p;

    pNorm=norm(p);
    qNorm=norm(q);

    l=l+(pNorm/qNorm)^2 * (pNorm-D)/D ;
    fprintf("it %i |p|=%g \t |q|=%g \t l=%g\n",iLoop,pNorm,qNorm,l)

end


end