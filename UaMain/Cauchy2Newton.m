
function [r2,rForce,rWork,DuC2N,DvC2N,DhC2N,DlC2N]=Cauchy2Newton(func,step,Variables,nNodes,xCauchy,xNewton)



C2N=xCauchy+step*(xNewton-xCauchy);


if Variables=="-uvhl-"

    DuC2N=C2N(1:nNodes) ;
    DvC2N=C2N(nNodes+1:2*nNodes);
    DhC2N=C2N(2*nNodes+1:3*nNodes);
    DlC2N=C2N(3*nNodes+1:end);



    gamma=1;
    [r2,~,~,rForce,rWork]=func(gamma,DuC2N,DvC2N,DhC2N,DlC2N) ;

elseif Variables=="-uvl-"

    DuC2N=C2N(1:nNodes) ;
    DvC2N=C2N(nNodes+1:2*nNodes);
    DhC2N=[] ; 
    DlC2N=C2N(2*nNodes+1:end);



    gamma=1;
    [r2,~,~,rForce,rWork]=func(gamma,DuC2N,DvC2N,DlC2N) ;



else

    error("case not found")


end

return


end