
function [L,b]=createLh(Nnodes,hfixednode,hfixedvalue,htiedA,htiedB)

%%

if isempty(hfixednode) && isempty(hfixedvalue) && isempty(htiedA) && isempty(htiedB)
    L=[];
    b=[];
    return
end



if numel(hfixednode) ~= numel(hfixedvalue)
    numel(hfixednode)
    numel(hfixedvalue)
    error(' # of fixed nodes must be equal to # of fixed values ');
end


if numel(htiedA) ~= numel(htiedB)
    numel(htiedA)
    numel(htiedB)
    error(' # of master nodes must be equal to # of slave nodes ');
end

Nconstrains=numel(hfixednode)+numel(htiedA);
nL=numel(hfixednode)+2*numel(htiedA); % number of non-zero elements in L

ia=zeros(nL,1) ; ib=zeros(nL,1) ; xval=zeros(nL,1) ;
b=zeros(Nconstrains,1);

II=1; JJ=1;

% h ties
for I=1:numel(htiedA) 
    ia(II)=JJ ; ib(II)=htiedA(I) ; xval(II)=1; II=II+1;
    ia(II)=JJ ; ib(II)=htiedB(I) ; xval(II)=-1; II=II+1;
    b(JJ)=0;  JJ=JJ+1;
end

% h fixed
for I=1:numel(hfixednode)
    ia(II)=JJ ; ib(II)=hfixednode(I) ; xval(II)=1 ; b(JJ)=hfixedvalue(I) ; II=II+1 ; JJ=JJ+1;
end

L=sparse(ia,ib,xval,Nconstrains,Nnodes);

end

