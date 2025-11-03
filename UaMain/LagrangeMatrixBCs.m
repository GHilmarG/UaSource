function [Luv,Luvrhs,lambdauv,Lh,Lhrhs,lambdah,CtrlVar]=...
        LagrangeMatrixBCs(CtrlVar,MUA,ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,hfixednode,hfixedvalue,htiedA,htiedB,FixedNormalVelocityNode,FixedNormalVelocityValue)
    
    
    
    if numel(utiedA) ~= numel(utiedB) ; save TestSave ; error(' number of elements  in utiedA and utiedB  not the same ') ; end
    if numel(vtiedA) ~= numel(vtiedB) ; save TestSave ; error(' number of elements  in utiedA and utiedB not the same ') ; end
    
    % [----- BC
    % get rid of duplicate boundary conditions and just ignore extra BCs
    [ufixednodet,itemp]=unique(ufixednode) ; ufixedvalue=ufixedvalue(itemp);
    [vfixednodet,itemp]=unique(vfixednode) ; vfixedvalue=vfixedvalue(itemp);
    
    if numel(ufixednode) ~= numel(ufixednodet)  ; disp(' Duplicate Dirichlet BCs for u') ; end
    if numel(vfixednode) ~= numel(vfixednodet)  ; disp(' Duplicate Dirichlet BCs for v') ; end
    
    ufixednode=ufixednodet; vfixednode=vfixednodet;
    
    % velocity
    [Luv,Luvrhs]=CreateLuv(MUA,ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,FixedNormalVelocityNode,FixedNormalVelocityValue);
    
    
    lambdauv=Luvrhs*0;
    
    [Lh,Lhrhs]=createLh(MUA.Nnodes,hfixednode,hfixedvalue,htiedA,htiedB);
    
    lambdah=Lhrhs*0;
    
    
    %]
    
    
    CtrlVar.NouvTies=0 ; CtrlVar.NohTies=0 ;
    if isempty(utiedA) && isempty(vtiedB) && isempty(FixedNormalVelocityNode) ; CtrlVar.NouvTies=1  ;end
    if isempty(htiedA) ; CtrlVar.NohTies=1 ; end
    
    
    if CtrlVar.NouvTies ;
        [m,n]=size(Luv); if ~isequal(Luv*Luv',sparse(1:m,1:m,1)) ; save TestSave ; error(' Luv transpose(Luv) expected to be equal to the identity matrix, but is not!') ; end
    end
    
    if CtrlVar.NohTies ;
        [m,n]=size(Lh); if ~isequal(Lh*Lh',sparse(1:m,1:m,1)) ; save TestSave  ; error(' Lh transpose(Lh) expected to be equal to the identity matrix, but is not!') ; end
    end
       
    Luv=CtrlVar.BCsWeights*Luv;
    Luvrhs=CtrlVar.BCsWeights*Luvrhs;
    

    Lh=CtrlVar.BCsWeights*Lh;
    Lhrhs=CtrlVar.BCsWeights*Lhrhs;


    
    
end