function [hfixednode,hfixedvalue,htiedA,htiedB]=FromLhTohBCs(Lh,Lhrhs)


[nConstraints,nNodes]=size(Lh);
% vel constraints
[i,j,k]=find(Lh);

htiedA=zeros(nConstraints,1)+NaN; htiedB=zeros(nConstraints,1)+NaN;
hfixednode=zeros(nConstraints,1)+NaN; hfixedvalue=zeros(nConstraints,1)+NaN;

nhTies=0; nhFixed=0;
%for l=N(:,1)'
for l=1:nConstraints
    
    I=find(i==l) ;
    nI=numel(I);
    if nI==2
        if full(Lh(i(I(1)),j(I(1)))+Lh(i(I(2)),j(I(2))))==0   % it this really a tie?
            nodeA=j(I(1)) ;
            nodeB=j(I(2)) ;
            
            nhTies=nhTies+1;
            htiedA(nhTies)=nodeA;
            htiedB(nhTies)=nodeB;
            
        else
            fprintf('FromLhTohBCs: Found a multi-linear constrain on nodes (%i,%i) that is not just a simple tie\n',nodeA,nodeB)
        end
    elseif nI==1
        node=j(I(1));
        if full(Lh(i(I(1)),j(I(1)))==1)
            
            nhFixed=nhFixed+1;
            hfixednode(nhFixed)=node;
            hfixedvalue(nhFixed)=Lhrhs(l);
            
        else
            fprintf('FromLhTohBCs: Found a multi-linear constrain for node (%i) not representing a fixed nodal value \n',node)
        end
    else
        fprintf('FromLhTohBCs: Found a multi-linear constrain that is neither a nodal tie nor a fixed nodal value\n')
    end
    
end

I=~isnan(htiedA) ; htiedA=htiedA(I);  htiedB=htiedB(I); htiedA=htiedA(:) ; htiedB=htiedB(:) ;
I=~isnan(hfixednode) ; hfixednode=hfixednode(I);  hfixedvalue=hfixedvalue(I); hfixednode=hfixednode(:) ; hfixedvalue=hfixedvalue(:) ;

if numel(htiedA)==0 ; htiedA=[]; end
if numel(htiedB)==0 ; htiedB=[]; end
if numel(hfixednode)==0 ; hfixednode=[];end
if numel(hfixedvalue)==0 ; hfixedvalue=[];end

end

