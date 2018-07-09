

function [Lambda,Mu,l]=SolveAdjointEquationsDiscrete(kv,Lambda,Mu,BoundaryNodes,l,coordinates,connectivity,nip,...
        hModel,sModel,uModel,vModel,uMeas,vMeas,wMeasInt,...
        Nnodes,icount,iMisfitType,M)
    
    % Note: have still to indtroduce the error weighting matrix M :   B'  M (B u-d)
    
%     disp(' Step#2: Solve the adjoint model (right-hand side is a vector) ')
    
    utiedA=[] ; utiedB=[];  vtiedA=[] ; vtiedB=[];
    vfixednode=BoundaryNodes;  vfixedvalue=BoundaryNodes*0;
    ufixednode=BoundaryNodes;  ufixedvalue=BoundaryNodes*0;
    [L,Lb]=CreateLuv(Nnodes,ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,[],[]);
    
    
    if iMisfitType==10
         bModel=sModel-hModel;
         [~,B] = VertVelMatrixVector(hModel,bModel,uModel,vModel,coordinates,connectivity,nip); % Note; Here I should be using BModel instead of bModel
         % M is  Nele*nip x 2*Nodes
         
         
    
         % M is 2*Nnodes+Nele*nip x 2*Nodes
         Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
         k=2*Nnodes ; l=Nele*nip;
         B=[sparse(1:k,1:k,1,k,k) ; B];
         
         %r=[uModel-uMeas;vModel-vMeas;wModelInt-wMeasInt] ;
         %rh=2*B'*r;
         
         rh=2*B'*M*(B*[uModel;vModel]-[uMeas;vMeas;wMeasInt]);
    else
        rh=2*[uModel-uMeas;vModel-vMeas] ;
    end
    
    [sol,l]=solveKApeSymmetric(kv',L,rh,Lb,[Lambda;Mu],l,icount);
    Lambda=sol(1:Nnodes) ; Mu=sol(Nnodes+1:2*Nnodes);
    
    
end
