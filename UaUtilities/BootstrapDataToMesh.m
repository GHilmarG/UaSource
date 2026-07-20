


function [MeasAtNodes,StandardDeviation]=BootstrapDataToMesh(CtrlVar,MUA,Meas,nStride,nBootstrapReplicates,options)

%%
%
% Simple bootstrapping estimate of interpolation errors when interpolating 2D spatial data set (measurements) to nodal locations.
%
% Inputs:
%
% Meas.x     : x location-of measurements
% Meas.y     : y location-of measurements
% Meas.Value : measurements
%
% nstride  : stride between measurements, i.e. if nstride=10, only each 10th measurement will be used.
% nBootstrapReplicates : number of bootstrap replicates.
%
%
% Outputs:
%
% MeasAtNodes       : interpolated and averaged bootstrap values at the nodes of the mesh
% StandardDeviation : standard deviation of bootstrap values at each nodal location.
%
%
% This works fine to come up with an rough estimation of interpolation errors as one interpolates scattered set of
% measurements onto the nodes of the mesh. But if number of replicates is low, estimated standard deviation at nodal
% locations where a measurement is already available within a short distance might easily be an underestimate.
%
% This approach does not properly deal with spatial correlation in data.
%
%
%%

arguments (Input)
    CtrlVar struct
    MUA     struct
    Meas {mustBeA(Meas,["struct","numeric"])}
    nStride {mustBeNumeric} = 1
    nBootstrapReplicates {mustBeNumeric} = 10
    options.Plot logical = false
end


%% First get rid of all measurements located outside of FE mesh boundary
fprintf("BringDataToMesh: Getting rid of measurements outside of mesh boundary.\n")
io=inpoly2([Meas.x Meas.y], [MUA.Boundary.x MUA.Boundary.y]);
Meas.x=Meas.x(io) ; Meas.y=Meas.y(io) ;  Meas.Value=Meas.Value(io);

Meas.x=Meas.x(1:nStride:end);
Meas.y=Meas.y(1:nStride:end);
Meas.Value=Meas.Value(1:nStride:end);

nMeas=numel(Meas.Value);

%% Now use simple bootstrapping, interpolating data onto nodes using bootstrap replicates to estimate uncertainty at query locations.

%nBootstrapReplicates=20 ;

BootStrapEstimates=nan(MUA.Nnodes,nBootstrapReplicates) ;
rng(0);

xNodes=MUA.coordinates(:,1);
yNodes=MUA.coordinates(:,2);

fprintf("Bootstrapping estimation:\n")

% add a bit of jitter, so that the scattered interpolant does not average repeated values.

if ~isfield(MUA,"EleAreas")
    MUA.EleAreas=TriAreaFE(MUA.coordinates,MUA.connectivity); % areas for each element
    MUA.Area=sum(MUA.EleAreas);                               % total FE mesh area
end


ds=sqrt(min(MUA.EleAreas)) ;
xjitter=0.01*ds*randn(nMeas,1);
yjitter=0.01*ds*randn(nMeas,1);

parfor I=1:nBootstrapReplicates

    idx=randi(nMeas,nMeas,1);


    xMeas=Meas.x(idx)+xjitter;
    yMeas=Meas.y(idx)+yjitter;
    ValMeas=Meas.Value(idx);
    F=scatteredInterpolant(xMeas,yMeas,ValMeas);
    BootStrapEstimates(:, I) = F(xNodes,yNodes);


end

MeasAtNodes = mean(BootStrapEstimates, 2, 'omitnan');
StandardDeviation = std(BootStrapEstimates, 0, 2, 'omitnan');


%%
if options.Plot
    Fig=FindOrCreateFigure("Bootstrap");
    T=tiledlayout("flow");

    T1=nexttile;
    cbar=UaPlots(CtrlVar,MUA,[],MeasAtNodes,CreateNewFigure=false);
    title("Mean of measurements")
    subtitle("")
    hold on
    %plot(Measurements.x/1000,Measurements.y/1000,".k")
    title(cbar,"(m)")


    T2=nexttile;
    cbar=UaPlots(CtrlVar,MUA,[],StandardDeviation,CreateNewFigure=false);
    set(gca,'ColorScale','lin')
    %clim([0.1 100])
    hold on;
    plot(Meas.x/1000,Meas.y/1000,".k",MarkerSize=10)
    plot(xNodes/1000,yNodes/1000,".r",MarkerSize=10)
    title("Standard deviation of measurements")
    subtitle("")
    title(cbar,"(m)")
    %%

    %% examples of some histogram of estimates at few selected nodes
    iNodes=randi(MUA.Nnodes,5,1) ;

    for k=1:5
        iNodes(k);
        FindOrCreateFigure("node"+num2str(k))
        histogram(BootStrapEstimates(iNodes(k),:));
        title(sprintf("node=%i    std=%g",iNodes(k),StandardDeviation(iNodes(k))))
    end




end
end






