
function [etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=...
    calcStrainRatesEtaInt(CtrlVar,MUA,ub,vb,AGlen,n)

% [etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,u,v,AGlen,n)
% calculates strain rates and effective viscosity at integration points and Eint needed for the directional
% derivative
%
% If AGlen is empty, the effective viscosity is not calculated and a NaN is returned
%
% vectorized
%
%[etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,ub,vb,AGlen,n);

narginchk(4,6)

if nargin<5   || isempty(AGlen)
    etaInt=[] ; Eint=[] ;  txx=[] ; tyy=[] ; txy=[] ;
end

if nargin>5
    if ~isempty(n)
        if CtrlVar.AGlenisElementBased
            n=n+zeros(MUA.Nele,1);
        else
            n=n+zeros(MUA.Nnodes,1);
        end
    end
end


ndim=2;


ubnod=reshape(ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vbnod=reshape(vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
coox=reshape(MUA.coordinates(MUA.connectivity,1),MUA.Nele,MUA.nod);
cooy=reshape(MUA.coordinates(MUA.connectivity,2),MUA.Nele,MUA.nod);

if nargin>4 && ~isempty(AGlen)
    if CtrlVar.AGlenisElementBased
        AGlennod=[];
        nGlennod=[];
    else
        AGlennod=reshape(AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
        nGlennod=reshape(n(MUA.connectivity,1),MUA.Nele,MUA.nod);
    end
end

exx=zeros(MUA.Nele,MUA.nip); eyy=zeros(MUA.Nele,MUA.nip); exy=zeros(MUA.Nele,MUA.nip);
xint=zeros(MUA.Nele,MUA.nip) ; yint=zeros(MUA.Nele,MUA.nip);

if nargin>4
    if ~isempty(AGlen)
        AGlenint=zeros(MUA.Nele,MUA.nip);
        nGlenint=zeros(MUA.Nele,MUA.nip);
    end
end



for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    
    xint(:,Iint)=coox*fun;
    yint(:,Iint)=cooy*fun;
    
    if nargin>4
        if ~isempty(AGlen)
            if CtrlVar.AGlenisElementBased
                AGlenint(:,Iint)=AGlen;
                nGlenint(:,Iint)=n;
            else
                temp=AGlennod*fun;
                temp(temp<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
                AGlenint(:,Iint)=temp;
                nGlenint(:,Iint)=nGlennod*fun;
            end
        end
    end
    
    if exist('MUA','var') && isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
        Deriv=MUA.Deriv(:,:,:,Iint);
    else
        Deriv=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,MUA.points,Iint);
    end
    % [Deriv]=derivVector(coordinates,connectivity,nip,Iint); % Nele x dof x nod
    % The derivative depends on the det of each element
    % is therefore a much bigger array than fun
    
    
    
    for I=1:MUA.nod
        exx(:,Iint)=exx(:,Iint)+Deriv(:,1,I).*ubnod(:,I);
        eyy(:,Iint)=eyy(:,Iint)+Deriv(:,2,I).*vbnod(:,I);
        exy(:,Iint)=exy(:,Iint)+0.5*(Deriv(:,1,I).*vbnod(:,I) + Deriv(:,2,I).*ubnod(:,I));
    end
end

e=real(sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2));

if nargin>4
    if ~isempty(AGlen)
        n=nGlenint;
        % etaInt=real(0.5*AGlenint.^(-1./n).*e.^((1-n)./n));
        % Eint=real((1-n)./(4*n).*AGlenint.^(-1./n).*e.^((1-3*n)./n));

        [etaInt,Eint,e,dEtadA]=EffectiveViscositySSTREAM(CtrlVar,AGlenint,n,exx,eyy,exy) ;
    end
end

if nargin>4 && nargout>8
    if ~isempty(AGlen)
        txx=2*etaInt.*exx ;
        tyy=2*etaInt.*eyy ;
        txy=2*etaInt.*exy ;
    end
end

if nargin>4
    if any(isnan(etaInt)) ; save TestSave  ;  error(' NaN in etaInt ' ) ; end
    if any(isnan(Eint)) ; save TestSave  ;  error(' NaN in Eint ' ) ; end
end

if any(isnan(e)) ; save TestSave  ;  error(' NaN in effective strain rates  ' ) ; end


% Directional derivative is
%
% Eint (   (2 exx+eyy) d Delta u/dx  + (2 eyy +exx) d Delta v/dy  + exy (d Delta v/ dx + d Delta u /dy))
%


end
