
%% Get some data

%load('D:\GoogleDriveStreamingHomeHP\My Drive\Runs\Calving\PIG-TWG\ResultsFiles\0214000-FR2020to2500-20km-uvh-Tri3-SlidWeertman-Duvh-MRIM6HadGEM2-abMask0A-P-BCVel-kH10000-TM0k2-Alim-Clim-Ca1-Cs100000-Aa1-As100000-VelITS120-BM3-SMB_RACHMO2k3_2km-.mat')

load('D:\GoogleDriveStreamingHomeHP\My Drive\Runs\Calving\PIG-TWG\ResultsFiles\0233640-FR2020to2500-10km-uvh-Tri3-SlidWeertman-Duvh-MRlASE3-abMask0A-IOR-P-BCVel-kH10000-TM0k2-Alim-Clim-Ca1-Cs100000-Aa1-As100000-VelITS120-BM3-SMB_RACHMO2k3_2km-.mat')

% load('D:\GoogleDriveStreamingHomeHP\My Drive\Runs\Calving\PIG-TWG\ResultsFiles\0210340-FR2020to2500-5km-uvh-Tri3-SlidWeertman-Duvh-MRIM6HadGEM2-abMask0A-P-BCVel-kH10000-TM0k2-Alim-Clim-Ca1-Cs100000-Aa1-As100000-VelITS120-BM3-SMB_RACHMO2k3_2km-.mat')

% load('D:\GoogleDriveStreamingHomeHP\My Drive\Runs\Calving\PIG-TWG\ResultsFiles\0205140-FR2020to2500-2k5km-uvh-Tri3-SlidWeertman-Duvh-MRIM6HadGEM2-abMask0A-P-BCVel-kH10000-TM0k2-Alim-Clim-Ca1-Cs100000-Aa1-As100000-VelITS120-BM3-SMB_RACHMO2k3_2km-.mat')

%
% 
% $$q=\int_{\partial \mathcal{A}} \rho_i h \, \mathbf{v} \cdot \mathbf{n} \, d \Gamma = \int_{\mathcal{A}}   \nabla \cdot (\rho_i h \mathbf{v} ) \, d \mathcal{A} $$
% 

%% line integral

tic
[qGL,qGLx,qGLy,Fub,Fvb,Fr,Fh,LakeNodes,GLgeo]=FluxAcrossGroundingLine(CtrlVar,MUA,F.GF,F.ub,F.vb,F.ud,F.vd,F.h,F.rho);

FluxLineIntegral=sum(qGL) ; 
toc


FindOrCreateFigure("GL Flux") ; 
PlotMuaMesh(CtrlVar,MUA) ; hold on ; 
plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'g','LineWidth',2); 
scale=1; hold on ;  quiver(GLgeo(:,7)/CtrlVar.PlotXYscale,GLgeo(:,8)/CtrlVar.PlotXYscale,qGLx,qGLy,scale,'color','r') ; 
axis equal


% Now do the area integral approach 
tic
qx=F.rho.*F.h.*F.ub;
qy=F.rho.*F.h.*F.vb;

[dqxdx,dqxdy]=calcFEderivativesMUA(qx,MUA,CtrlVar) ; 
[dqydx,dqydy]=calcFEderivativesMUA(qy,MUA,CtrlVar) ; 

Dint=dqxdx+dqydy ;

Dnod=ProjectFintOntoNodes(MUA,Dint);

Dnod=Dnod.*F.GF.node;

Int=FEintegrate2D(CtrlVar,MUA,Dnod); 

FluxAreaIntegral=sum(Int) ;
toc


E=100*(abs(FluxAreaIntegral-FluxLineIntegral)) /FluxLineIntegral ; % relative difference in %

fprintf("Ground line flux: \n")
fprintf("\t  Line integral estimate: %g (10^9 kg/yr) \n",FluxLineIntegral/1e9)
fprintf("\t  Area integral estimate: %g (10^9 kg/yr) \n",FluxAreaIntegral/1e9)
fprintf("\t Abs Relative Difference: %g%% \n",E)

%%