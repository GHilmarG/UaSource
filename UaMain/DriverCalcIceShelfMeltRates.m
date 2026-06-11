function ab=DriverCalcIceShelfMeltRates(MUA,VelDataSet,nSmooth,minab,maxab)
    
    %%
    % VelDataSet: '450m'|'990m'
    
    CtrlVar=Ua2D_DefaultParameters();
    time=0;
    [xGL,yGL]=ReadBindschadlerGroundingLine;
    
    %CtrlVar.ReadInitialMeshFileName='MeshFileEle118037Nod6.mat';  % PIG and Thwaites
    
    % in Runs/AntIceShelfThinning
    %CtrlVar.ReadInitialMeshFileName='MeshFileEle18722Node6';  % whole of Ant
    %CtrlVar.ReadInitialMeshFileName='MeshFileEle363209Nod6';  % whole of Ant fine
    
    %MUA=CreateMUA(CtrlVar,connectivity,coordinates);
    
    %CtrlVar.ReadInitialMeshFileName='MeshFileEle118037Nod6.mat';
    
    %load(CtrlVar.ReadInitialMeshFileName) ;
    
    if nargin<2
        VelDataSet='450m';
    end
    
    if nargin<3 || isempty(nSmooth)
        nSmooth=0; 
    end
    
    if nargin< 4
        minab=[] ; maxab=[];
    end
    
    load('Prognostic-RestartFile','ub','vb','MUA')
    x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
%    [u,v,Err,Mask]=EricVelocities(CtrlVar,MUA.coordinates,VelDataSet);
    u=ub ; v=vb; 
    
    figure ; QuiverColorGHG(x,y,u,v,CtrlVar);
    
    Experiment=[];
    CtrlVar.fidlog=1;
    [UserVar,s,b,S,B,alpha]=GetGeometry(UserVar,CtrlVar,MUA,time);

    
    h=s-b; time =0 ;
    
    [UserVar,rho,rhow,g]=GetDensities(UserVar,CtrlVar,MUA,time,s,b,h,S,B);
    GF=GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);
    GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
    [UserVar,as,ab]=GetMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
    dsdt=s*0;
   
    
    [Fdhdt,dhdt]=GetMeasuredIceShelfThinningRates(x,y);
    
    ab=CalcIceShelfMeltRates(CtrlVar,MUA,u,v,s,b,S,B,rho,rhow,dsdt,as,dhdt);
    
    
    
    if ~isempty(maxab)
        ab(ab>maxab)=maxab;
    end
    
    if ~isempty(minab)
        ab(ab<minab)=minab;
    end
    
    % smoothing by projecting onto elements and back to nodes
    [M,ElePerNode] = Ele2Nodes(MUA.connectivity,MUA.Nnodes);
    
    for Ismooth=1:nSmooth
        ab=Nodes2EleMean(MUA.connectivity,ab);  
        ab=M*ab;
    end
    
    dhdt(GF.node>0.5)=NaN;
    
    CtrlVar.PlotXYscale=1000;
    figure ; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,dhdt,CtrlVar); 
    title('dhdt')
    ylabel(colorbar,'dhdt (m/a)')
    
    
    ab(GF.node>0.5)=NaN;
    
    %%
    CtrlVar.PlotXYscale=1000;
    figure ; 
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,-ab,CtrlVar); 
    hold on 
    %plot(xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,'k');
    plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'k','LineWidth',1);
    title('Ice Shelf balance melt rate') 
    ylabel(colorbar,'melt rate (m/a)')
    xlabel('xps (km)')
    ylabel('yps (km)')
    %%
 
    
 %   BalanceMeltRate=ab ; connectivityBalanceMeltRate=connectivity ; coordinatesBalanceMeltRate=coordinates;   
 %   save('BalanceIceShelfMeltRates.mat','connectivityBalanceMeltRate','coordinatesBalanceMeltRate','BalanceMeltRate')
end
