function [Results,Info]=DhdtMisfitVersusBlockSize(CtrlVar,MUA,F,Ls,x0,y0,xPix,yPix,dPix,PixelSize,IncludeFunction,options)

%%
%
%   [Results,Info]=DhdtMisfitVersusBlockSize(CtrlVar,MUA,F,Ls,x0,y0,xPix,yPix,dPix,PixelSize,IncludeFunction)
%   [Results,Info]=DhdtMisfitVersusBlockSize(...,Err=Err,Plot=true)
%
% Diagnostic for the "-blocks-" dh/dt misfit (CtrlVar.Inverse.DataMisfit.dhdt.Method="-blocks-"): statistics of the residuals
% between modelled and measured dh/dt averaged over square blocks, as a function of the block size L.
%
% The modelled dh/dt is calculated once at the integration points (dhdtAtPoints.m), using the fields in F (u, v, h, as, ab, rho), so
% F must contain a forward (uv) solution, for example the final state of an inversion. For each block size, the blocks and the
% block-averaging operator are defined exactly as in CreateDhdtBlockMeasurements.m, with the same selection criteria, and the model
% and the measurements are averaged over the same blocks.
%
% The results show the scale above which the model reproduces the measurements, and how the residuals decrease with L. This can be
% used to select the block size L, and the errors of the block averages. If the residuals are dominated by errors that are
% uncorrelated at scales above the smallest L, the RMS of the residuals decreases as 1/L (Info.Slope=-1).
%
% Inputs: as for CreateDhdtBlockMeasurements.m, except
%
%   Ls              : vector of block sizes, each a multiple of PixelSize (L=PixelSize gives the pixel footprints)
%
% Optional (name=value):
%
%   Err   : errors of the block averages, used to calculate chi2=2J/N. Either a scalar, or a function handle of the form Err(L),
%           Err(xc,yc) or Err(xc,yc,L), where (xc,yc) are the block centres. Default: [] (chi2 not calculated).
%   Plot  : if true, the RMS of the residuals, and the standard deviations of the measured and the modelled block averages, are
%           plotted as functions of L. Default: false.
%
% Outputs:
%
%   Results : table with one row per block size and the columns
%
%       L             block size
%       nBlocks       number of blocks retained
%       MeanResidual  mean of the residuals (model-measurement)
%       RMSResidual   RMS of the residuals
%       StdObs        standard deviation of the measured block averages
%       StdModel      standard deviation of the modelled block averages
%       Chi2          2J/N=mean((residual/Err)^2), NaN if Err is not given. A value of about one indicates that Err is consistent
%                     with the residuals.
%
%   Info    : structure with Info.Slope, the slope of log(RMSResidual) against log(L) (least-squares fit over the block sizes with
%             at least 3 blocks), and Info.Blocks{i}, the block centres, measured and modelled block averages for the i-th block size.
%
% *Example*
%
%   Ls=[10 20 40 60 80 100]*1e3 ;
%   x0=gx(1)-PixelSize/2 ; y0=gy(1)-PixelSize/2 ;         % block edges on pixel edges
%   Grounded=@(x,y) MaskAntarctica(x,y)==2 ;
%   [Results,Info]=DhdtMisfitVersusBlockSize(CtrlVar,MUA,F,Ls,x0,y0,xPix,yPix,dPix,PixelSize,Grounded,Err=1,Plot=true)
%
% see also: CreateDhdtBlockMeasurements.m, dhdtObservationMisfit.m, dhdtAtPoints.m
%
%%

arguments
    CtrlVar struct
    MUA struct
    F
    Ls (:,1) double
    x0 (1,1) double
    y0 (1,1) double
    xPix double
    yPix double
    dPix double
    PixelSize (1,1) double
    IncludeFunction
    options.Err = []
    options.Plot (1,1) logical = false
end

% modelled dh/dt at the integration points, calculated once
Pint=ShapeFunctionsAtIntegrationPoints(CtrlVar,MUA);
hdot=dhdtAtPoints(MUA,F,Pint);

CtrlVarQuiet=CtrlVar;
CtrlVarQuiet.InfoLevel=0;

nL=numel(Ls);
nBlocks=zeros(nL,1) ; MeanResidual=nan(nL,1) ; RMSResidual=nan(nL,1) ; StdObs=nan(nL,1) ; StdModel=nan(nL,1) ; Chi2=nan(nL,1) ;
Info.Blocks=cell(nL,1);

for i=1:nL

    L=Ls(i);
    M=CreateDhdtBlockMeasurements(CtrlVarQuiet,MUA,struct(),L,x0,y0,xPix,yPix,dPix,PixelSize,IncludeFunction,1);
    nBlocks(i)=numel(M.dhdtData);

    if nBlocks(i)==0
        continue
    end

    Model=M.dhdtA*hdot;
    Obs=M.dhdtData;
    r=Model-Obs;

    MeanResidual(i)=mean(r);
    RMSResidual(i)=sqrt(mean(r.^2));
    StdObs(i)=std(Obs);
    StdModel(i)=std(Model);

    if ~isempty(options.Err)
        Err=options.Err;
        if isa(Err,'function_handle')
            switch nargin(Err)
                case 1
                    s=Err(L);
                case 2
                    s=Err(M.dhdtx,M.dhdty);
                otherwise
                    s=Err(M.dhdtx,M.dhdty,L);
            end
        else
            s=Err;
        end
        s=s(:)+zeros(nBlocks(i),1);
        Chi2(i)=mean((r./s).^2);
    end

    Info.Blocks{i}=struct('x',M.dhdtx,'y',M.dhdty,'Obs',Obs,'Model',Model);

end

Results=table(Ls,nBlocks,MeanResidual,RMSResidual,StdObs,StdModel,Chi2, ...
    'VariableNames',{'L','nBlocks','MeanResidual','RMSResidual','StdObs','StdModel','Chi2'});

I=nBlocks>=3 & isfinite(RMSResidual) & RMSResidual>0;
if sum(I)>=2
    pfit=polyfit(log(Ls(I)),log(RMSResidual(I)),1);
    Info.Slope=pfit(1);
else
    Info.Slope=NaN;
end

if options.Plot
    fig=FindOrCreateFigure("dh/dt misfit versus block size") ; clf(fig) ;
    loglog(Ls/1000,RMSResidual,'o-',Ls/1000,StdObs,'s-',Ls/1000,StdModel,'d-') ;
    hold on
    if any(I)
        iRef=find(I,1);
        loglog(Ls/1000,RMSResidual(iRef)*Ls(iRef)./Ls,'k--') ;   % slope -1 for reference
    end
    xlabel("Block size, $L$ (km)",Interpreter="latex") ;
    ylabel("$\dot{h}$ (m/yr)",Interpreter="latex") ;
    legend("RMS(model$-$measured)","std(measured)","std(model)","$\propto 1/L$",Interpreter="latex",Location="best") ;
    title("Block averages of $\dot{h}$, residual statistics versus block size",Interpreter="latex") ;
    subtitle(sprintf("slope of log(RMS) against log(L): %.2f",Info.Slope)) ;
    grid on
end

end
