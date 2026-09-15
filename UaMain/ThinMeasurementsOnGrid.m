

function [Thinned,Info]=ThinMeasurementsOnGrid(x,y,d,ds,Options)

%%
%
%   [Thinned,Info]=ThinMeasurementsOnGrid(x,y,d,ds,Options)
%
% Thins densely sampled scattered measurements by snapping them onto a regular grid of
% spacing ds and averaging within each occupied cell. Returns a reduced data set,
% together with an error estimate for each thinned measurement based on the variability
% of the raw data within its cell.
%
% This is mesh free. The output is an ordinary data set, with locations and values, and
% can be used with any mesh, at any resolution, and saved and reused indefinitely. That
% is the main reason for preferring it to binning against the elements of a particular
% FE mesh.
%
% Only occupied cells are created. Empty cells never appear and cost nothing, so ds can
% be made small without penalty on a sparse survey. For data collected along profiles,
% where the along-track sampling is far denser than the across-track line spacing, the
% reduction in the number of measurements is typically large.
%
%
%% Why thin at all
%
% Closely spaced measurements along a profile are not independent. Treating them as
% though they were makes a densely sampled profile count as many more independent
% constraints than it really provides, which over-weights the data against the prior
% and invites the inversion to fit structure that is not resolvable. Thinning to at or
% above the error correlation length removes that redundancy at essentially no cost,
% and avoids having to introduce a correlated error covariance at the measurement
% points, which would be dense and expensive.
%
%
%% Choosing ds
%
% Two requirements, pulling in opposite directions:
%
%   ds  should be at or above the correlation length of the measurement error, so that
%       the thinned measurements are genuinely close to independent, and
%
%   ds  should be below the smallest element size of any mesh the data will be used
%       with, so that the thinning, and not the mesh, never becomes the factor limiting
%       resolution.
%
% If the error correlation length exceeds the element size, these cannot both be met,
% and thinning alone will not do the job. A correlated error covariance would then be
% needed. It is worth checking the empirical semivariogram of the data before settling
% on ds, if only to confirm that this is not the case.
%
%
%% Placement of the thinned measurement
%
% The thinned measurement is placed at the centroid of the raw measurements within the
% cell, and NOT at the centre of the cell.
%
% Averaging the data over a cell and assigning the mean to a single location is exact
% only if the underlying field is linear across the cell. Placing the point at the
% centroid removes the first-order error exactly, leaving only a curvature term of
% order ds^2. Placing it at the cell centre would leave a first-order error wherever
% the measurements are distributed asymmetrically within the cell, which for data along
% profiles is essentially always.
%
% The residual second-order error is uniform in space, being set by ds, rather than
% varying with element size as it would for mesh-based binning. It therefore cannot
% be mistaken for structured misfit.
%
%
%% The error estimate
%
%   sigma^2 = RawError^2/n + SigmaCorrelated^2 + sRough^2
%
% The first term is the reduction of uncorrelated error on averaging n samples. The
% second is a prescribed floor that does not average down, which is where any error
% common to a whole profile, such as a navigation or datum offset, belongs. Without it
% a cell holding very many samples would be assigned an implausibly small error.
%
% sRough is a representativeness term: the within-cell scatter about the local linear
% trend, in excess of that already explained by RawError. It captures bed structure at
% scales below ds that a single thinned measurement cannot represent.
%
% The local trend is removed before the scatter is measured, so that a genuine and
% perfectly resolvable slope across the cell is not misread as noise.
%
%
%% Shrinkage, and why it is needed here
%
% With ds below the element scale most cells hold only a handful of measurements, and
% RSS/dof from three or four samples is an extremely noisy estimate of a variance. Used
% directly it would produce erratic sigma with no physical basis.
%
% The local estimate is therefore shrunk toward a pooled estimate formed over all
% cells,
%
%   s2 = (RSS + ShrinkageDof*s2Pooled)/(dof + ShrinkageDof)
%
% which is the posterior mean under an inverse-gamma prior with ShrinkageDof degrees of
% freedom. Cells with many samples are barely affected; cells with few are pulled
% toward the global value, which is what one would want in the absence of local
% evidence. Set Options.ShrinkageDof=0 to disable.
%
%
%% Rank of the local trend fit
%
% Measurements along a single profile crossing a cell are collinear, and a plane fit is
% then rank deficient across track. The rank is tested per cell and the fit degrades:
%
%   n>=4, well conditioned : plane fit,             dof = n-3
%   n>=3, collinear        : line fit along the
%                            principal direction,   dof = n-2
%   n==2                   : scatter about mean,    dof = 1
%   n==1                   : no local estimate, the pooled value is used
%
%
%% Inputs
%
%   x, y, d  : m x 1 raw measurement locations and values. NaN entries are discarded.
%
%   ds       : grid spacing, in the units of x and y.
%
%   Options  : optional structure
%
%     .RawError         []     Error of an individual raw measurement, scalar or m x 1.
%                              If empty, it is estimated from the pooled within-cell
%                              scatter and returned in Info.RawError. Note that such an
%                              estimate is an upper bound, since it also contains any
%                              real bed structure at scales below ds.
%
%     .SigmaCorrelated  0      Error floor that does not average down.
%
%     .ShrinkageDof     4      Strength of the shrinkage toward the pooled variance.
%
%     .SigmaFloor       0      Lower bound applied to the returned sigma.
%
%     .Origin           [0 0]  Grid origin. Fix this if several data sets are to be
%                              thinned onto a common grid.
%
%     .InfoLevel        1      Set to 0 to suppress the summary.
%
%
%% Outputs
%
%   Thinned  : structure with, for each occupied cell,
%                .x,.y        centroid of the raw measurements in the cell
%                .d           weighted mean of the measurements
%                .sigma       error estimate
%                .Count       number of raw measurements
%                .Scatter     sqrt of the shrunken within-cell variance
%                .Roughness   sRough, the excess over RawError
%                .FitType     3=plane, 2=line, 1=mean, 0=single sample
%                .xGrid,.yGrid  the snapped grid location, for reference
%                .ds, .Origin
%
%   Info     : structure with the pooled estimates and a summary of the reduction.
%
%
%% Example
%
%   [Thinned,Info]=ThinMeasurementsOnGrid(Bobs.x,Bobs.y,Bobs.B,250) ;
%   save("ThinnedBedPicks.mat","Thinned","Info")
%
% and later, on any mesh,
%
%   [O,Inside]=BuildNode2DataMap(CtrlVar,MUA,Thinned.x,Thinned.y) ;
%   r=O(Inside,:)*F.B-Thinned.d(Inside) ;
%   J=sum((r./Thinned.sigma(Inside)).^2)/2 ;
%
%
%% See also
%
% BuildNode2DataMap, BinMeasurementsOntoMesh
%
%%

narginchk(4,5)

if nargin<5 || isempty(Options)
    Options=struct ;
end

Default.RawError=[] ;
Default.SigmaCorrelated=0 ;
Default.ShrinkageDof=4 ;
Default.SigmaFloor=0 ;
Default.Origin=[0 0] ;
Default.InfoLevel=1 ;

Fields=fieldnames(Default) ;
for k=1:numel(Fields)
    if ~isfield(Options,Fields{k})
        Options.(Fields{k})=Default.(Fields{k}) ;
    end
end

x=x(:) ; y=y(:) ; d=d(:) ;
m=numel(x) ;

if numel(y)~=m || numel(d)~=m
    error("ThinMeasurementsOnGrid:InconsistentInputs", ...
        "x, y and d must all have the same number of elements.")
end

if ~isscalar(ds) || ~(ds>0)
    error("ThinMeasurementsOnGrid:BadSpacing","ds must be a positive scalar.")
end

RawError=Options.RawError ;
EstimateRawError=isempty(RawError) ;

if ~EstimateRawError
    RawError=RawError(:) ;
    if isscalar(RawError)
        RawError=RawError*ones(m,1) ;
    elseif numel(RawError)~=m
        error("ThinMeasurementsOnGrid:InconsistentInputs", ...
            "Options.RawError must be empty, a scalar, or have as many elements as x.")
    end
    if any(RawError<=0)
        error("ThinMeasurementsOnGrid:NonPositiveError","Options.RawError must be positive.")
    end
end

%% 1) Discard non-finite data

Keep=isfinite(x) & isfinite(y) & isfinite(d) ;

x=x(Keep) ; y=y(Keep) ; d=d(Keep) ;
if ~EstimateRawError
    RawError=RawError(Keep) ;
end

n0=numel(x) ;

if n0==0
    error("ThinMeasurementsOnGrid:NoData","No finite measurements were supplied.")
end

%% 2) Snap onto the grid, and identify the occupied cells
%
% Only occupied cells are created, so a fine grid over a sparse survey costs nothing.

x0=Options.Origin(1) ; y0=Options.Origin(2) ;

ix=round((x-x0)/ds) ;
iy=round((y-y0)/ds) ;

[C,~,bin]=unique([ix iy],"rows") ;

nBins=size(C,1) ;

Count=accumarray(bin,1,[nBins 1]) ;

%% 3) Cell centroids and weighted means

if EstimateRawError
    w=ones(n0,1) ;              % provisional, uniform weighting
else
    w=1./RawError.^2 ;
end

wSum=accumarray(bin,w,[nBins 1]) ;

xm=accumarray(bin,w.*x,[nBins 1])./wSum ;
ym=accumarray(bin,w.*y,[nBins 1])./wSum ;
dm=accumarray(bin,w.*d,[nBins 1])./wSum ;

%% 4) Within-cell scatter about the local linear trend
%
% Centred on the weighted centroid, so that the constant term drops out.

dx=x-xm(bin) ; dy=y-ym(bin) ; dd=d-dm(bin) ;

Sxx=accumarray(bin,dx.*dx,[nBins 1]) ;
Sxy=accumarray(bin,dx.*dy,[nBins 1]) ;
Syy=accumarray(bin,dy.*dy,[nBins 1]) ;
Sxd=accumarray(bin,dx.*dd,[nBins 1]) ;
Syd=accumarray(bin,dy.*dd,[nBins 1]) ;
Sdd=accumarray(bin,dd.*dd,[nBins 1]) ;

Det=Sxx.*Syy-Sxy.^2 ;
Tr =Sxx+Syy ;

WellConditioned = Det > 1e-8*Tr.^2 & Tr>0 ;

RSS=Sdd ;
dof=Count-1 ;
FitType=ones(nBins,1) ;
FitType(Count==1)=0 ;

% (a) plane fit
Plane = WellConditioned & Count>=4 ;

if any(Plane)
    b=( Syy(Plane).*Sxd(Plane)-Sxy(Plane).*Syd(Plane))./Det(Plane) ;
    c=(-Sxy(Plane).*Sxd(Plane)+Sxx(Plane).*Syd(Plane))./Det(Plane) ;
    RSS(Plane)=Sdd(Plane)-b.*Sxd(Plane)-c.*Syd(Plane) ;
    dof(Plane)=Count(Plane)-3 ;
    FitType(Plane)=3 ;
end

% (b) line fit along the principal direction, the usual case for data along a profile
Line = ~Plane & Count>=3 & Tr>0 ;

if any(Line)

    Lmax=(Tr+sqrt(max(Tr.^2-4*Det,0)))/2 ;

    v1=Sxy ; v2=Lmax-Sxx ;
    Swap=abs(Sxy)<=eps*Tr ;
    v1(Swap)=double(Sxx(Swap)>=Syy(Swap)) ;
    v2(Swap)=double(Sxx(Swap)< Syy(Swap)) ;

    vn=hypot(v1,v2) ; vn(vn==0)=1 ;
    v1=v1./vn ; v2=v2./vn ;

    u=dx.*v1(bin)+dy.*v2(bin) ;

    Suu=accumarray(bin,u.*u,[nBins 1]) ;
    Sud=accumarray(bin,u.*dd,[nBins 1]) ;

    Ok=Line & Suu>0 ;

    RSS(Ok)=Sdd(Ok)-Sud(Ok).^2./Suu(Ok) ;
    dof(Ok)=Count(Ok)-2 ;
    FitType(Ok)=2 ;

end

RSS=max(RSS,0) ;
dof=max(dof,0) ;

%% 5) Pooled variance, and shrinkage of the local estimates
%
% RSS/dof from three or four samples is a very poor estimate of a variance. The local
% estimates are therefore shrunk toward the pooled value.

TotalDof=sum(dof) ;

if TotalDof>0
    s2Pooled=sum(RSS)/TotalDof ;
else
    s2Pooled=0 ;
end

kShrink=Options.ShrinkageDof ;

s2=(RSS+kShrink*s2Pooled)./max(dof+kShrink,realmin) ;

if kShrink==0
    s2=zeros(nBins,1) ;
    Has=dof>=1 ;
    s2(Has)=RSS(Has)./dof(Has) ;
    s2(~Has)=s2Pooled ;
end

%% 6) The raw measurement error
%
% If not prescribed, the pooled within-cell scatter is used. This is an upper bound on
% the true pick error, since it also contains any real bed structure at scales below
% ds, and the representativeness term is then zero for all but the rougher cells.

if EstimateRawError
    RawVar=s2Pooled*ones(nBins,1) ;
    RawErrorOut=sqrt(s2Pooled) ;
else
    RawVar=accumarray(bin,RawError.^2,[nBins 1])./Count ;
    RawErrorOut=Options.RawError ;
end

%% 7) The error estimate

if EstimateRawError
    RawMeanVar=RawVar./Count ;          % uniform weights
else
    RawMeanVar=1./wSum ;                % inverse-variance weighted mean
end

sRough2=max(s2-RawVar,0) ;

sigma=sqrt(RawMeanVar+Options.SigmaCorrelated^2+sRough2) ;
sigma=max(sigma,Options.SigmaFloor) ;

%% 8) Output

Thinned.x=xm ;
Thinned.y=ym ;
Thinned.d=dm ;
Thinned.sigma=sigma ;
Thinned.Count=Count ;
Thinned.Scatter=sqrt(s2) ;
Thinned.Roughness=sqrt(sRough2) ;
Thinned.FitType=FitType ;
Thinned.xGrid=x0+C(:,1)*ds ;
Thinned.yGrid=y0+C(:,2)*ds ;
Thinned.ds=ds ;
Thinned.Origin=Options.Origin ;

Info.nRaw=m ;
Info.nFinite=n0 ;
Info.nThinned=nBins ;
Info.ReductionFactor=n0/nBins ;
Info.RawError=RawErrorOut ;
Info.RawErrorEstimated=EstimateRawError ;
Info.PooledScatter=sqrt(s2Pooled) ;
Info.PooledDof=TotalDof ;
Info.ds=ds ;

if Options.InfoLevel>=1
    fprintf("ThinMeasurementsOnGrid: ds=%g,  %i raw measurements -> %i thinned (factor %.1f). \n", ...
        ds,n0,nBins,Info.ReductionFactor)
    fprintf("                        samples per cell : min=%i, median=%g, max=%i \n", ...
        min(Count),median(Count),max(Count))
    if EstimateRawError
        fprintf("                        raw error estimated from pooled scatter : %g  (%i dof) \n", ...
            RawErrorOut,TotalDof)
    end
    fprintf("                        sigma            : min=%g, median=%g, max=%g \n", ...
        min(sigma),median(sigma),max(sigma))
    fprintf("                        fit types        : plane=%i, line=%i, mean=%i, single=%i \n", ...
        sum(FitType==3),sum(FitType==2),sum(FitType==1),sum(FitType==0))
end

end
