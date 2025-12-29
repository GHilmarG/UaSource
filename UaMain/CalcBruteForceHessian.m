



function [Hsparse,Hfull,g0,J0] = CalcBruteForceHessian(func,p0,CtrlVar,iRange)



%%
%
% p0 : are the model parameters (A,B,C)
%
% func : A function handle
%
%
% [J,g,H]=func(p)
%
% returns the cost (J), gradient (G) and Hessian (H). (Only cost, J, is used here.)
%
% The function returns a finite-differences estimate of dJ/dp for selected elements of p in the range iRange
%
% On return dJ is a vector of
%
% [dJ/dA(iRange) dJ/dB(iRange) dJ/dC(iRange)]
%
% assuming one is inverting for A, B and C.
%
%%



fprintf(' Calculating Hessian using brute-force method...')

if nargin < 4 || isempty(iRange) || anynan(iRange)
    iRange=1:100 ;
end

[J0,g0]=func(p0);

deltaStep=CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize*abs(p0);

if any(deltaStep==0)
    fprintf("At least one of the deltaStep values is equal to zero. \n")
    error("CalcBruteForceHessian:ZeroDelta","At least one of the deltaStep values is equal to zero.")
end


n=numel(p0);

Hfull=nan(n,n);

nColumns=numel(iRange) ;
Hvalues=nan(n*nColumns,1) ;


Loop="-parfor-";

HessianFiniteDifferenceType="central second order" ;
HessianFiniteDifferenceType="central fourth order" ;

if contains(Loop,"-parfor-")  % the parfor option can only be used when calculating ALL columns of the Hessian
    parfevalOnAll(@warning,0,'off','all');

    % I found it difficult to do the parfor over subset of columns, so here it is done for all of them
    tStart=tic;
    parfor k=1:n

        Hcol=ParForBody(k,func,p0,deltaStep(k),HessianFiniteDifferenceType) ;
        Hfull(:,k)=Hcol;
        

    end
    tEnd=toc(tStart);
    fprintf("parfor evaluation of the Hessian in %f sec",tEnd)
end



if contains(Loop,"-for-")
tic
for k=1:nColumns

    iColumn=iRange(k);

    p1=p0;
    p1(iColumn)=p1(iColumn)+deltaStep(iColumn); % change one element within p

    pm1=p0;
    pm1(iColumn)=pm1(iColumn)-deltaStep(iColumn);

    [J1,g1]=func(p1);
    [Jm1,gm1]=func(pm1);


    Hcol=(g1-gm1)/(2*deltaStep(iColumn)) ;

    Hfull(:,iColumn)=Hcol ;
    
    % Here I need to introduce some sparsity idea
    ind=(k-1)*n+1:(k-1)*n+n ;
    Hvalues(ind)=Hcol;

     if rem(k-1,nColumns/10)==0
         fprintf("%i ",round(100*k/nColumns))
     end

end

toc
end


iH=repmat(1:n,1,numel(iRange)) ; iH=iH(:) ;
jH=repmat(iRange,n,1) ; jH=jH(:);
Hsparse=sparse(iH,jH,Hvalues,n,n) ;

Hfull=0.5*(Hfull+Hfull');

fprintf("\n")




%  k=3 ; ind=((((k-1)*n+1):(k*n))') ; [iH(ind) jH(ind) Hvalues(ind)]
%  I=abs(H)>1; figure(100) ; hold off ; spy(I)
%
% FindOrCreateFigure("Hfull") ; contourf(rot90(abs(Hfull))) ; colorbar ;   set(gca,'ColorScale','log')
%
%  FindOrCreateFigure("Hfull sparsity") ; H=abs(Hfull) > mean(abs(Hfull)) ; spy(H)
% norm(Hfull-full(Hsparse))
%
% H \dp = -g 
%
end


function Hcol=ParForBody(k,func,p0,delta,HessianFiniteDifferenceType)


switch HessianFiniteDifferenceType

    case "central second order"

        %% central second order

        p1=p0;
        p1(k)=p1(k)+delta; % change one element within p

        pm1=p0;
        pm1(k)=pm1(k)-delta;

        [J1,g1]=func(p1);
        [Jm1,gm1]=func(pm1);


        Hcol=(g1-gm1)/(2*delta) ;

        %Hfull(:,iColumn)=Hcol ;

        % Here I need to introduce some sparsity idea
        %ind=(k-1)*n+1:(k-1)*n+n ;
        %Hvalues(ind)=Hcol;

    case "central fourth order"
        %% central fourth order
        p1=p0;
        pm1=p0;
        p2=p0;
        pm2=p0;
        p1(k)=p1(k)+delta;
        p2(k)=p2(k)+2*delta;
        pm1(k)=pm1(k)-delta;
        pm2(k)=pm2(k)-2*delta;

        [J1,g1]=func(p1);
        [J2,g2]=func(p2);
        [Jm1,gm1]=func(pm1);
        [Jm2,gm2]=func(pm2);

        Hcol=(gm2/12-2*gm1/3+2*g1/3-g2/12)/delta;

    otherwise

        error("case not found ")


end
end