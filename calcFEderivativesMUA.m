function [dfdx,dfdy,xint,yint]=calcFEderivativesMUA(f,MUA,CtrlVar)

% [dfdx,dfdy,xint,yint]=calcFEderivatives(f,MUA,CtrlVar)
% calculates x and y derivatives of a nodal variable at integration points



ndim=2;
[points,weights]=sample('triangle',MUA.nip,ndim);

fnod=reshape(f(MUA.connectivity,1),MUA.Nele,MUA.nod);

dfdx=zeros(MUA.Nele,MUA.nip); dfdy=zeros(MUA.Nele,MUA.nip);


% f is a vector with nod values
% the derivative at a given integration point is
% dfds=Dx f   ( [Nele x nod] * [nod]
% Dx=Deriv(:,1,:)  which is Nele x nod
% dfdx(nEle)=Dx

for Iint=1:MUA.nip
    
    %    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
    Deriv=MUA.Deriv(:,:,:,Iint);
    %    else
    %        Deriv=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    %    end

    
    for I=1:MUA.nod
        dfdx(:,Iint)=dfdx(:,Iint)+Deriv(:,1,I).*fnod(:,I);
        dfdy(:,Iint)=dfdy(:,Iint)+Deriv(:,2,I).*fnod(:,I);
    end
    
end

if nargout>2
    
    xint=zeros(MUA.Nele,MUA.nip) ; yint=zeros(MUA.Nele,MUA.nip);
    coox=reshape(MUA.coordinates(MUA.connectivity,1),MUA.Nele,nod);
    cooy=reshape(MUA.coordinates(MUA.connectivity,2),MUA.Nele,nod);
    
    for Iint=1:MUA.nip
        fun=shape_fun(Iint,ndim,nod,points) ;
        
        xint(:,Iint)=coox*fun;
        yint(:,Iint)=cooy*fun;
        
    end
else
    xint=[] ; yint=[];
end


end
