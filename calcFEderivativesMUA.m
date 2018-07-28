function [dfdx,dfdy,xint,yint]=calcFEderivativesMUA(f,MUA,CtrlVar)

%%
%
%   [dfdx,dfdy,xint,yint]=calcFEderivativesMUA(f,MUA,CtrlVar)
% 
% calculates x and y derivatives of a nodal variable at integration points
%
%   f          :  a nodal quantity 
%   dfdx, dfdy :  x and y derivatitves of f 
%   xint, yint :  x, y locations of the elements of dfdx and dfdy 
%
%   CtrlVar    : optional, can be left empty. 
% 
% Note: On input f is a nodal variable, ie defined at nodes
%       On return dfdx and dfdy are integration-point variables, ie defined at
%       the integration points xint and yint.
%
% See also: ProjectFintOntoNodes
%%


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
    
    Deriv=MUA.Deriv(:,:,:,Iint);
        
    for I=1:MUA.nod
        dfdx(:,Iint)=dfdx(:,Iint)+Deriv(:,1,I).*fnod(:,I);
        dfdy(:,Iint)=dfdy(:,Iint)+Deriv(:,2,I).*fnod(:,I);
    end
    
end

if nargout>2
    
    xint=zeros(MUA.Nele,MUA.nip) ; yint=zeros(MUA.Nele,MUA.nip);
    coox=reshape(MUA.coordinates(MUA.connectivity,1),MUA.Nele,MUA.nod);
    cooy=reshape(MUA.coordinates(MUA.connectivity,2),MUA.Nele,MUA.nod);
    
    for Iint=1:MUA.nip
        fun=shape_fun(Iint,ndim,MUA.nod,points) ;
        
        xint(:,Iint)=coox*fun;
        yint(:,Iint)=cooy*fun;
        
    end
else
    xint=[] ; yint=[];
end


end
