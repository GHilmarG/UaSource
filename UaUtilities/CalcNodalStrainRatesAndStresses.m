function [tbx,tby,txx,tyy,txy,exx,eyy,exy,e,eta]=CalcNodalStrainRatesAndStresses(CtrlVar,UserVar,MUA,F)

narginchk(4,4)

%%
% Calculates strains and devitoric stresses.
%
% [txzb,tyzb,txx,tyy,txy,exx,eyy,exy,e,eta]=CalcNodalStrainRatesAndStresses(CtrlVar,UserVar,MUA,F) 
%
% Strains and stresses are first calculated at integration points, then projeted onto nodes.
%
% On output, all variables are nodal variables.
%
%   txzb, tyzb    : x and y components of the basal shear stresses (i.e. not x and y components of basal traction).
%   txx,tyy,txy   : horizontal deviatoric stresses
%   exx,eyy,exy   : horizontal strain rates
%   e             : effective strain rate
%   eta           : effective viscosity
%
% the basal stress caculation is done using the basal boundary condition as:
%
%   txzb = tbx + ( 2 txx + tyy) \p_x b + txy \p_y b
%       < N_p | N_q >  txzb_q = < N_p | tbx + ( 2 txx + tyy) \p_x b + txy \p_y b >
%
% Cauchy stresses can then be calculated as \sigma_{xx}=2 \tau_{xx} + \tau_{yy} + \sigma_{zz}
% where \sigma_{zz}= - \rho g (s-z)
%
% Upper surface stresses are \sigma_{xx}=2 \tau_{xx} + \tau_{yy}
% Lower surface stresses are \sigma_{xx}=2 \tau_{xx} + \tau_{yy} - \rho g h
%
%
% Example:
% 
%   load('CrackRestartfileExample.mat','CtrlVarInRestartFile','MUA','F','BCs','GF')
%   CtrlVar=CtrlVarInRestartFile;
%   UserVar=UserVarInRestartFile
%   [txzb,tyzb,txx,tyy,txy,exx,eyy,exy,e,eta]=CalcNodalStrainRatesAndStresses(CtrlVar,UserVar,MUA,F);
%   x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2); 
%   [X,Y]=ndgrid(linspace(min(x),max(x),20),linspace(min(y),max(y),20));
%   I=nearestNeighbor(MUA.TR,[X(:) Y(:)]);  % find nodes within computational grid closest to the regularly scape X and Y grid points.
%   scale=1e-3;
%   FigStrainAndStresses=figure; 
%   PlotTensor(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,txx(I),txy(I),tyy(I),scale);
%   hold on
%   PlotMuaBoundary(CtrlVar,MUA,'k')
%   axis equal
%
%
%%


[tbx,tby,tb] = CalcBasalTraction(CtrlVar,UserVar,MUA,F); % returns nodal values and uses the sliding law 
[etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,F.ub,F.vb,F.AGlen,F.n); % returns integration point values

ndim=2; neq=MUA.Nnodes;

bnod=reshape(F.b(MUA.connectivity,1),MUA.Nele,MUA.nod);
snod=reshape(F.s(MUA.connectivity,1),MUA.Nele,MUA.nod);


% 
% % tbxnod=reshape(tbx(MUA.connectivity,1),MUA.Nele,MUA.nod);
% % tbynod=reshape(tby(MUA.connectivity,1),MUA.Nele,MUA.nod);
% 
% % [points,weights]=sample('triangle',MUA.nip,ndim);
% 
% 
% Tx=zeros(MUA.Nele,MUA.nod);
% Ty=zeros(MUA.Nele,MUA.nod);
% 
% 
% % vector over all elements for each integration point
% for Iint=1:MUA.nip
% 
%     fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
% 
%     if ~(isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ))
%         fprintf("CalcNodalStrainRatesAndStresses: MUA is not in an expected state. MUA updated. Consider doing this ahead of this call. \n")
%         fprintf("             MUA=UpdateMUA(CtrlVar,MUA)   \n")
%         MUA=UpdateMUA(CtrlVar,MUA) ;
%     end
% 
%     Deriv=MUA.Deriv(:,:,:,Iint);
%     detJ=MUA.DetJ(:,Iint);
% 
% 
%     dsdx=zeros(MUA.Nele,1); dsdy=zeros(MUA.Nele,1);
%     dbdx=zeros(MUA.Nele,1); dbdy=zeros(MUA.Nele,1);
% 
%     % derivatives for all elements at this integration point
%     for Inod=1:MUA.nod
%         dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
%         dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
%         dbdx=dbdx+Deriv(:,1,Inod).*bnod(:,Inod);
%         dbdy=dbdy+Deriv(:,2,Inod).*bnod(:,Inod);
%     end
% 
%     dbdx=kk_proj(dbdx,-CtrlVar.dbdxZero,CtrlVar.dbdxZero);
%     dbdy=kk_proj(dbdy,-CtrlVar.dbdyZero,CtrlVar.dbdyZero);
% 
% 
%     tbxint=tbxnod*fun;  % values at this integration point
%     tbyint=tbynod*fun;
% 
%     txzint=tbxint+(2*txx(:,Iint)+tyy(:,Iint)).*dbdx+txy(:,Iint).*dbdy;
%     tyzint=tbyint+txy(:,Iint).*dbdx+(2*tyy(:,Iint)+txx(:,Iint)).*dbdy;
% 
%     detJw=detJ*MUA.weights(Iint);
% 
%     for Inod=1:MUA.nod
% 
%         Tx(:,Inod)=Tx(:,Inod)+txzint.*fun(Inod).*detJw;
%         Ty(:,Inod)=Ty(:,Inod)+tyzint.*fun(Inod).*detJw;
% 
% 
%     end
% end
% 
% % assemble right-hand side
% 
% rhx=sparseUA(neq,1); rhy=sparseUA(neq,1);
% for Inod=1:MUA.nod
%     rhx=rhx+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Tx(:,Inod),neq,1);
%     rhy=rhy+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Ty(:,Inod),neq,1);
% end
% M=MassMatrix2D1dof(MUA);
% sol=M\[rhx rhy] ;
% txzb=full(sol(:,1)) ; tyzb=full(sol(:,2));

if nargout>2
    [txx,tyy,txy,exx,eyy,exy,e,eta]=ProjectFintOntoNodes(MUA,txx,tyy,txy,exx,eyy,exy,e,etaInt);
end

% if ~isreal(txzb) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:txzbNotReal','txzb not real!') ; end
% if ~isreal(tyzb) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:tyzbNotReal','tyzb not real!') ; end
if ~isreal(txx) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:txxbNotReal','txx not real!') ; end
if ~isreal(txx) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:txxbNotReal','txx not real!') ; end
if ~isreal(tyy) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:tyybNotReal','tyy not real!') ; end
if ~isreal(txy) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:tyybNotReal','txy not real!') ; end
if ~isreal(exx) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:exxbNotReal','exx not real!') ; end
if ~isreal(eyy) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:eyybNotReal','eyy not real!') ; end
if ~isreal(exy) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:eyybNotReal','exy not real!') ; end
if ~isreal(e) ; save TestSave ; error('CalcNodalStrainRatesAndStresses:ebNotReal','e not real!') ; end



end


