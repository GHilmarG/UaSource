

function GenerateFieldWithGivenCovarianceInTwoDimentions(x,y,sigma,ampl)

%%
%
% Creates a 2D field with a given covariance. Just a sketch, wrote this quickly, can easily be optimized further...
%
%
% 
%
%%

if nargin==0
    % generate some example data
    nPoints=7000;
    xy=rand(nPoints,2);
    x=xy(:,1) ;
    y=xy(:,2);

    nugget=0.1;
    sigma=0.05; %horizontal correlation length

    ampl=5;

end

dist=zeros(numel(x),numel(y));

CovarianceMatrix=zeros(numel(x),numel(y)) ;

for ix=1:numel(x)
    for jx=ix:numel(y)
        dist(ix,jx)= sqrt((x(ix)-x(jx))^2 + (y(ix)-y(jx))^2);
        dist(jx,ix)=dist(ix,jx) ;
    end
end


for ix=1:numel(x)
    for jx=ix:numel(y)
        CovarianceMatrix(ix,jx)=ampl*exp(-dist(ix,jx)/sigma) ;
        CovarianceMatrix(jx,ix)=CovarianceMatrix(ix,jx);
    end
end


L=chol(CovarianceMatrix,"lower");


ErrorNormal=randn(numel(x),1); 

ErrorCorrelated=L*ErrorNormal;

ErrorCorrelated=ErrorCorrelated+sqrt(nugget)*ErrorNormal;

%% Plot
DT=delaunayTriangulation(x,y);

figure(1000) ;
hold off
PlotNodalBasedQuantities(DT.ConnectivityList,DT.Points,ErrorCorrelated) ;
hold on
plot(x,y,".r")
axis equal tight
CM=cmocean('-balanced',25,'pivot',0) ; colormap(CM);


%% semivariogram 

dd=sigma/50;
nd=0.5/dd ;

[sv,sv_dist,nav]=get_sv(x,y,ErrorCorrelated,dd,nd) ;

figure(2000) ; 
plot(sv_dist,sv,"o-",LineWidth=2)
title("Semivariogram ")
xlabel("distance")
ylabel("sv")
xline(sigma,"--","sigma")   
% I think the amplitude estimate might assuming that the nugget effect is small
%yline(2*(ampl+nugget),"--","2*(amplitude+nugget)","LabelHorizontalAlignment","left")   
yline(ampl+nugget,"--","amplitude+nugget","LabelHorizontalAlignment","left")   

yline(nugget,"--","nugget")   
yL=ylim; ylim([0 yL(2)])
%%



end