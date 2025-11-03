
function [Luv,b]=CreateLuv(MUA,ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,FixedNormalVelocityNode,FixedNormalVelocityValue)

if isempty(ufixednode) && isempty(ufixedvalue) && isempty(vfixednode) && isempty(vfixedvalue) ...
        && isempty(utiedA) && isempty(utiedB) && isempty(vtiedA) && isempty(vtiedB) ...
        && isempty(FixedNormalVelocityNode)
    Luv=[];
    b=[];
    return
end

if any(isnan(ufixedvalue)) ; error(' nan in ufixedvalue ') ; end
if any(isnan(vfixedvalue)) ; error(' nan in vfixedvalue ') ; end

if numel(utiedA) ~= numel(utiedB)
    numel(utiedA)
    numel(utiedB)
    error(' # of master nodes must be equal to # of slave nodes ');
end

if numel(vtiedA) ~= numel(vtiedB)
    numel(vtiedA)
    numel(vtiedB)
    error(' # of master nodes must be equal to # of slave nodes ');
end


Nconstrains=numel(ufixednode)+numel(vfixednode)+numel(utiedA)+numel(vtiedA)+numel(FixedNormalVelocityNode);
nL=numel(ufixednode)+numel(vfixednode)+2*numel(utiedA)+2*numel(vtiedA)+2*numel(FixedNormalVelocityNode); % number of non-zero elements in L


ia=zeros(nL,1) ; ib=zeros(nL,1) ; xval=zeros(nL,1) ;
b=zeros(Nconstrains,1);

II=1; JJ=1;
% JJ is the constraint number, ie the equation number in Luv (each equation is represented by one row in the matrix)

% u ties
for I=1:numel(utiedA) 
    ia(II)=JJ ; ib(II)=utiedA(I) ; xval(II)=1; II=II+1;
    ia(II)=JJ ; ib(II)=utiedB(I) ; xval(II)=-1; II=II+1;
    b(JJ)=0;  JJ=JJ+1;
end

% v ties
for I=1:numel(vtiedA) 
    
    ia(II)=JJ ; ib(II)=vtiedA(I)+MUA.Nnodes ; xval(II)=1; II=II+1;
    ia(II)=JJ ; ib(II)=vtiedB(I)+MUA.Nnodes ; xval(II)=-1; II=II+1;
    b(JJ)=0; JJ=JJ+1;
end


% u fixed
for I=1:numel(ufixednode)
    ia(II)=JJ ; ib(II)=ufixednode(I) ; xval(II)=1 ; b(JJ)=ufixedvalue(I) ; II=II+1 ; JJ=JJ+1;
end


% v fixed
for I=1:numel(vfixednode)
    ia(II)=JJ ; ib(II)=vfixednode(I)+MUA.Nnodes ; xval(II)=1 ; b(JJ)=vfixedvalue(I) ; II=II+1 ; JJ=JJ+1;
end

if numel(FixedNormalVelocityNode)>0
    
    [nx,ny,xn,yn,Nx,Ny] = CalcEdgeAndNodalNormals(MUA.connectivity,MUA.coordinates,MUA.Boundary.Edges);
    for I=1:numel(FixedNormalVelocityNode)
        ia(II)=JJ ; ib(II)=FixedNormalVelocityNode(I)         ; xval(II)=Nx(FixedNormalVelocityNode(I)); II=II+1;
        ia(II)=JJ ; ib(II)=FixedNormalVelocityNode(I)+MUA.Nnodes  ; xval(II)=Ny(FixedNormalVelocityNode(I)); II=II+1;
        b(JJ)=-FixedNormalVelocityValue(I);
        JJ=JJ+1;
    end
end

Luv=sparseUA(ia,ib,xval,Nconstrains,2*MUA.Nnodes);

end

