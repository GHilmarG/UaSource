function [txzb,tyzb,txx,tyy,txy,exx,eyy,exy,e,eta]=CalcNodalStrainRatesAndStresses(CtrlVar,MUA,AGlen,n,C,m,GF,s,b,ub,vb,ud,vd)

% Calculates strains and devitoric stresses. 
%
% [txzb,tyzb,txx,tyy,txy,exx,eyy,exy,e,eta]=CalcNodalStrainRatesAndStresses(CtrlVar,MUA,AGlen,n,C,m,GF,s,b,ub,vb,ud,vd)
%
% Strains and stresses are first calculated at integration points, then projeted onto nodes.
% On output all variables are nodal variables. 
%
% txzb, tyzb    : x and y components of the basal shear stresses (i.e. not x and y components of basal traction).
% txx,tyy,txy   : horizontal deviatoric stresses 
% exx,eyy,exy   : horizontal strain rates
% e             : effective strain rate  
% eta           : effective viscosity
%
% the basal stress caculation is done using the basal boundary condition as:
% txzb = tbx + ( 2 txx + tyy) \p_x b + txy \p_y b
%   < N_p | N_q >  txzb_q = < N_p | tbx + ( 2 txx + tyy) \p_x b + txy \p_y b >
%
% Stresses can then be calculated as \sigma_{xx}=2 \tau_{xx} + \tau_{yy} + \sigma_{zz}
% wheere \sigma_{zz}= - \rho g (s-z)
% Upper surface stresses are \sigma_{xx}=2 \tau_{xx} + \tau_{yy} 
% Lower surface stresses are \sigma_{xx}=2 \tau_{xx} + \tau_{yy} - \rho g h 
%
%

%narginchk(11,13)
%nargoutchk(1,9);



%save CalcNodalStrainRatesAndStressesTestSave

%%

[tbx,tby,tb,beta2] = CalcBasalTraction(CtrlVar,MUA,ub,vb,C,m,GF); % returns nodal values
[etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,ub,vb,AGlen,n); % returns integration point values

ndim=2; neq=MUA.Nnodes;

bnod=reshape(b(MUA.connectivity,1),MUA.Nele,MUA.nod);
snod=reshape(s(MUA.connectivity,1),MUA.Nele,MUA.nod);



tbxnod=reshape(tbx(MUA.connectivity,1),MUA.Nele,MUA.nod);
tbynod=reshape(tby(MUA.connectivity,1),MUA.Nele,MUA.nod);

[points,weights]=sample('triangle',MUA.nip,ndim);


Tx=zeros(MUA.Nele,MUA.nod);
Ty=zeros(MUA.Nele,MUA.nod);


% vector over all elements for each integartion point
for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ')
        Deriv=MUA.Deriv(:,:,:,Iint);
        detJ=MUA.DetJ(:,Iint);
    else
        [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end
        
    dsdx=zeros(MUA.Nele,1); dsdy=zeros(MUA.Nele,1);
    dbdx=zeros(MUA.Nele,1); dbdy=zeros(MUA.Nele,1);

    % derivatives for all elements at this integration point
    for Inod=1:MUA.nod
        dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
        dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
        dbdx=dbdx+Deriv(:,1,Inod).*bnod(:,Inod);
        dbdy=dbdy+Deriv(:,2,Inod).*bnod(:,Inod);
    end
    
    dbdx=kk_proj(dbdx,-CtrlVar.dbdxZero,CtrlVar.dbdxZero);
    dbdy=kk_proj(dbdy,-CtrlVar.dbdyZero,CtrlVar.dbdyZero);

    
    tbxint=tbxnod*fun;  % values at this integration point
    tbyint=tbynod*fun;
    
    txzint=tbxint+(2*txx(:,Iint)+tyy(:,Iint)).*dbdx+txy(:,Iint).*dbdy;
    tyzint=tbyint+txy(:,Iint).*dbdx+(2*tyy(:,Iint)+txx(:,Iint)).*dbdy;
     
   detJw=detJ*weights(Iint);
    
    for Inod=1:MUA.nod
        
        Tx(:,Inod)=Tx(:,Inod)+txzint.*fun(Inod).*detJw;
        Ty(:,Inod)=Ty(:,Inod)+tyzint.*fun(Inod).*detJw;
        
        
    end
end

% assemble right-hand side

rhx=sparse2(neq,1); rhy=sparse2(neq,1);
for Inod=1:MUA.nod
    rhx=rhx+sparse2(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Tx(:,Inod),neq,1);
    rhy=rhy+sparse2(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Ty(:,Inod),neq,1);
end


M=MassMatrix2D1dof(MUA);

sol=M\[rhx rhy] ; 
txzb=full(sol(:,1)) ; tyzb=full(sol(:,2));

if nargout>2
    [txx,tyy,txy,exx,eyy,exy,e,eta]=ProjectFintOntoNodes(MUA,txx,tyy,txy,exx,eyy,exy,e,etaInt);
end

if ~isreal(txzb) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:txzbNotReal','txzb not real!') ; end
if ~isreal(tyzb) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:tyzbNotReal','tyzb not real!') ; end
if ~isreal(txx) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:txxbNotReal','txx not real!') ; end
if ~isreal(txx) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:txxbNotReal','txx not real!') ; end
if ~isreal(tyy) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:tyybNotReal','tyy not real!') ; end
if ~isreal(txy) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:tyybNotReal','txy not real!') ; end
if ~isreal(exx) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:exxbNotReal','exx not real!') ; end
if ~isreal(eyy) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:eyybNotReal','eyy not real!') ; end
if ~isreal(exy) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:eyybNotReal','exy not real!') ; end
if ~isreal(e) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:ebNotReal','e not real!') ; end
    


end


