function [dfdx,dfdy,xint,yint]=calcFEderivatives(f,coordinates,connectivity,nip,CtrlVar)
  
    % [dfdx,dfdy,xint,yint]=calcFEderivatives(f,coordinates,connectivity,nip,CtrlVar)
    % calculates x and y derivatives of a nodal variable at integration points
    % There is a better routine for doing this [dfdx,dfdy,xint,yint]=calcFEderivativesMUA(f,MUA,CtrlVar)
    
    
    
    [Nele,nod]=size(connectivity); ndim=2;
    [points,weights]=sample('triangle',nip,ndim);
    
    fnod=reshape(f(connectivity,1),Nele,nod);
    
    coox=reshape(coordinates(connectivity,1),Nele,nod);
    cooy=reshape(coordinates(connectivity,2),Nele,nod);
    
    dfdx=zeros(Nele,nip); dfdy=zeros(Nele,nip); 
    xint=zeros(Nele,nip) ; yint=zeros(Nele,nip); 

    % f is a vector with nod values
    % the derivative at a given integration point is
    % dfds=Dx f   ( [Nele x nod] * [nod]
    % Dx=Deriv(:,1,:)  which is Nele x nod
    % dfdx(nEle)=Dx
    
    for Iint=1:nip                        
        fun=shape_fun(Iint,ndim,nod,points) ; 
        
        xint(:,Iint)=coox*fun;
        yint(:,Iint)=cooy*fun;
        [Deriv]=derivVector(coordinates,connectivity,nip,Iint); %  Deriv : Nele x dof x nod
        
        for I=1:nod
            dfdx(:,Iint)=dfdx(:,Iint)+Deriv(:,1,I).*fnod(:,I);
            dfdy(:,Iint)=dfdy(:,Iint)+Deriv(:,2,I).*fnod(:,I);
        end
        
    end
    
    
end
