function [ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,FixedNormalVelocityNode,FixedNormalVelocityValue]=...
    FromLuvTouvBCs(Luv,Luvrhs)

[nConstraints,nNodes]=size(Luv);
nNodes=nNodes/2; 
% vel constraints
[i,j,k]=find(Luv);



utiedA=zeros(nConstraints,1)+NaN; utiedB=zeros(nConstraints,1)+NaN;
vtiedA=zeros(nConstraints,1)+NaN; vtiedB=zeros(nConstraints,1)+NaN;
ufixednode=zeros(nConstraints,1)+NaN; vfixednode=zeros(nConstraints,1)+NaN;
ufixedvalue=zeros(nConstraints,1)+NaN; vfixedvalue=zeros(nConstraints,1)+NaN;
FixedNormalVelocityNode=zeros(nConstraints,1)+NaN; FixedNormalVelocityValue=zeros(nConstraints,1)+NaN;

nuTies=0; nvTies=0; nuFixed=0; nvFixed=0; nuvTies=0;
%for l=N(:,1)'
for l=1:nConstraints
    
    I=find(i==l) ;
    nI=numel(I);  % two non-zero entries in line l, ie two nodes involved in this constraint
    if nI==2 
        nodeA=j(I(1)) ;
        nodeB=j(I(2)) ;
        if full(Luv(i(I(1)),j(I(1)))+Luv(i(I(2)),j(I(2))))==0   % it this really a tie?
            
            if nodeA<=nNodes && nodeB<=nNodes
                nuTies=nuTies+1;
                utiedA(nuTies)=nodeA;
                utiedB(nuTies)=nodeB;
            elseif nodeA>nNodes && nodeB>nNodes
                nvTies=nvTies+1;
                vtiedA(nvTies)=nodeA-nNodes;
                vtiedB(nvTies)=nodeB-nNodes;
            end
        else
            if nodeA>nNodes ; nodeA=nodeA-nNodes ; end
            if nodeB>nNodes ; nodeB=nodeB-nNodes ; end
            
            if nodeA==nodeB   % u and v constraint on same node
                nuvTies=nuvTies+1;
                FixedNormalVelocityNode(nuvTies)=nodeA;
                FixedNormalVelocityValue(nuvTies)=Luvrhs(l);
            else
                fprintf('FromLuvTouvBCs: Found a multi-linear constrain on nodes (%i,%i) that is not just a simple tie\n',nodeA,nodeB)
            end
            
        end
    elseif nI==1
        node=j(I(1));
        if full(abs(Luv(i(I(1)),j(I(1))))==1)
            if node<=nNodes
                nuFixed=nuFixed+1;
                ufixednode(nuFixed)=node;
                ufixedvalue(nuFixed)=Luvrhs(l);
            elseif node>nNodes
                nvFixed=nvFixed+1;
                vfixednode(nvFixed)=node-nNodes;
                vfixedvalue(nvFixed)=Luvrhs(l);
            end
        else
            fprintf('FromLuvTouvBCs: Found a multi-linear constrain for node (%i) not representing a fixed nodal value \n',node)
        end
    else
        fprintf('FromLuvTouvBCs: Found a multi-linear constrain that is neither a nodal tie nor a fixed nodal value\n')
    end
    
end

I=~isnan(utiedA) ; utiedA=utiedA(I);  utiedB=utiedB(I); utiedA=utiedA(:) ; utiedB=utiedB(:) ;
I=~isnan(vtiedA) ; vtiedA=vtiedA(I);  vtiedB=vtiedB(I); vtiedA=vtiedA(:) ; vtiedB=vtiedB(:) ;

I=~isnan(ufixednode) ; ufixednode=ufixednode(I);  ufixedvalue=ufixedvalue(I); ufixednode=ufixednode(:) ; ufixedvalue=ufixedvalue(:) ;
I=~isnan(vfixednode) ; vfixednode=vfixednode(I);  vfixedvalue=vfixedvalue(I); vfixednode=vfixednode(:) ; vfixedvalue=vfixedvalue(:) ;

I=~isnan(FixedNormalVelocityNode) ; 
FixedNormalVelocityNode=FixedNormalVelocityNode(I);  FixedNormalVelocityValue=FixedNormalVelocityValue(I); 
FixedNormalVelocityNode=FixedNormalVelocityNode(:);  FixedNormalVelocityValue=FixedNormalVelocityValue(:) ;

if numel(utiedA)==0 ; utiedA=[]; end
if numel(utiedB)==0 ; utiedB=[]; end
if numel(vtiedA)==0 ; vtiedA=[]; end
if numel(vtiedA)==0 ; vtiedA=[]; end

if numel(ufixednode)==0 ; ufixednode=[]; end
if numel(vfixednode)==0 ; vfixednode=[]; end
if numel(ufixedvalue)==0 ; ufixedvalue=[]; end
if numel(vfixedvalue)==0 ; vfixedvalue=[]; end
if numel(FixedNormalVelocityNode)==0 ; FixedNormalVelocityNode=[]; end
if numel(FixedNormalVelocityValue)==0 ; FixedNormalVelocityValue=[]; end


end

