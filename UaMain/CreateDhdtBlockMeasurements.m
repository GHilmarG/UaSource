function [Meas,Info]=CreateDhdtBlockMeasurements(CtrlVar,MUA,Meas,L,x0,y0,xPix,yPix,dPix,PixelSize,IncludeFunction,Err)

%%
%
%   [Meas,Info]=CreateDhdtBlockMeasurements(CtrlVar,MUA,Meas,L,x0,y0,xPix,yPix,dPix,PixelSize,IncludeFunction,Err)
%
% Defines block-averaged measurements of dh/dt for the "-blocks-" dh/dt misfit,
%
%   CtrlVar.Inverse.DataMisfit.dhdt.Method="-blocks-" ;
%
% in which the modelled dh/dt is averaged over square blocks of size L, and compared with the average of a gridded dh/dt product over
% the same blocks (see dhdtObservationMisfit.m):
%
% $$ \bar{\dot{h}}_k = \frac{\sum_{p \in k} w_p \dot{h}_p}{\sum_{p \in k} w_p} , \qquad
%    J = \frac{1}{2} \sum_k \frac{(\bar{\dot{h}}_k - d_k)^2}{\sigma_k^2} $$
%
% where p are the integration points within block k, w_p their quadrature weights, and d_k the mean of the product pixels in block k.
% Each integration point is assigned to the block containing it, which approximates the block footprint to within the element size.
%
% Inputs:
%
%   L               : block size. Should be a multiple of the pixel size of the product.
%   x0, y0          : origin of the block grid, i.e. the block edges are at x0+i*L and y0+j*L. Should be on pixel edges, so that
%                     each block contains whole pixels.
%   xPix,yPix,dPix  : pixel centres and values of the gridded dh/dt product (NaN for missing values)
%   PixelSize       : pixel size of the product
%   IncludeFunction : function handle, Include=IncludeFunction(x,y), returning true where measurements are to be used, for example
%                     over grounded ice according to a fixed mask. Evaluated at the integration points. Use [] to include everything.
%                     The selection should not depend on quantities that change during the inversion.
%   Err             : errors (standard deviations) of the block averages, either a scalar or a function handle Err(xc,yc) of the
%                     block centres.
%
% A block is retained if (1) it lies fully within the mesh (tested at 5x5 points within the block), (2) at least 99% of its area is
% included according to IncludeFunction, and (3) all its product pixels have finite values.
%
% On return, Meas.dhdtx, Meas.dhdty (block centres), Meas.dhdtData, Meas.dhdtErr, Meas.dhdtA (block-averaging operator, acting on
% values at the integration points) and Meas.dhdtBlockSize are defined. Meas.dhdtA depends on the mesh and on the integration rule,
% and this function must be called again if either changes.
%
% Errors are treated as uncorrelated. As the misfit is likely to be dominated by model errors at short scales, L and Err should
% reflect the scale at which the model can be expected to reproduce the data.
%
% see also: CreateDhdtPointMeasurements.m, dhdtObservationMisfit.m, ShapeFunctionsAtIntegrationPoints.m
%
%%

narginchk(12,12)

MinIncluded=0.99;

xPix=xPix(:) ; yPix=yPix(:) ; dPix=dPix(:) ;

if abs(L/PixelSize-round(L/PixelSize))>1e-6
    warning('CreateDhdtBlockMeasurements:BlockSize','The block size L=%g is not a multiple of the pixel size %g.',L,PixelSize)
end
Offset=mod([xPix-x0 ; yPix-y0]-PixelSize/2,PixelSize);
if any(min(Offset,PixelSize-Offset)>1e-6*PixelSize)
    warning('CreateDhdtBlockMeasurements:Alignment','The block edges (x0,y0) are not aligned with the pixel edges.')
end

%% integration points and their blocks
Pint=ShapeFunctionsAtIntegrationPoints(CtrlVar,MUA);
bxy=[floor((Pint.x-x0)/L) floor((Pint.y-y0)/L)];
[Blocks,~,kb]=unique(bxy,'rows');
nb=size(Blocks,1);

Covered=accumarray(kb,Pint.w,[nb 1]);
if isempty(IncludeFunction)
    Include=true(size(Pint.x));
else
    Include=logical(IncludeFunction(Pint.x,Pint.y));
    Include=Include(:);
end
IncludedArea=accumarray(kb,Pint.w.*Include,[nb 1]);

%% (1) blocks fully within the mesh: test 5x5 points within each block
s=([0:4]'+0.5)/5;   %#ok<NBRAK1>
[sx,sy]=ndgrid(s,s);
xt=x0+(Blocks(:,1)+sx(:)')*L;  yt=y0+(Blocks(:,2)+sy(:)')*L;     % nb x 25
[~,InsideT]=BuildNode2DataMap(CtrlVar,MUA,xt(:),yt(:));
FullyInside=all(reshape(InsideT,nb,[]),2);

%% (3) product pixels in each block
pb=[floor((xPix-x0)/L) floor((yPix-y0)/L)];
[tf,loc]=ismember(pb,Blocks,'rows');
ok=tf & isfinite(dPix);
nPix=accumarray(loc(ok),1,[nb 1]);
SumPix=accumarray(loc(ok),dPix(ok),[nb 1]);
nPixExpected=round((L/PixelSize)^2);

Keep=FullyInside & IncludedArea>=MinIncluded*Covered & nPix==nPixExpected;
iKeep=find(Keep);
nKeep=numel(iKeep);

%% the block-averaging operator
NewIndex=zeros(nb,1); NewIndex(iKeep)=1:nKeep;
InKept=Keep(kb);
ip=find(InKept);
A=sparse(NewIndex(kb(ip)),ip,Pint.w(ip)./Covered(kb(ip)),nKeep,numel(Pint.w));

xc=x0+(Blocks(iKeep,1)+0.5)*L;
yc=y0+(Blocks(iKeep,2)+0.5)*L;

if isa(Err,'function_handle')
    BlockErr=Err(xc,yc);
else
    BlockErr=Err+zeros(nKeep,1);
end

Meas.dhdtx=xc;
Meas.dhdty=yc;
Meas.dhdtData=SumPix(iKeep)./nPix(iKeep);
Meas.dhdtErr=BlockErr(:);
Meas.dhdtA=A;
Meas.dhdtP=[];
Meas.dhdtBlockSize=L;

Info.nBlocksTouched=nb;
Info.nNotFullyInsideMesh=sum(~FullyInside);
Info.nNotIncluded=sum(FullyInside & IncludedArea<MinIncluded*Covered);
Info.nMissingPixels=sum(FullyInside & IncludedArea>=MinIncluded*Covered & nPix~=nPixExpected);
Info.nRetained=nKeep;
Info.CoveredAreaRatio=Covered(iKeep)/L^2;

if isempty(CtrlVar) || ~isfield(CtrlVar,'InfoLevel') || CtrlVar.InfoLevel>=1
    fprintf('CreateDhdtBlockMeasurements: L=%g, %i blocks retained (of %i; %i not fully inside the mesh, %i not included, %i with missing pixels).\n', ...
        L,Info.nRetained,nb,Info.nNotFullyInsideMesh,Info.nNotIncluded,Info.nMissingPixels)
end

if nKeep==0
    warning('CreateDhdtBlockMeasurements:NoData','No dh/dt blocks retained.')
end

end
