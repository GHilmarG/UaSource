
function [EleError,NormalizedEleError,rmsEleError]=fe2dErrorEstimate(iRefinement,Experiment,coordinates0,connectivity0,...
        Boundary0,MeshBoundaryCoordinates,u0,v0,h0,s0,b0,AGlen0,C0,n,m,rho,rhow,g,Itime,nip,DTxy0,TRIxy0,...
        CtrlVar)
        

% estimates errors for each element


x0=coordinates0(:,1) ; y0=coordinates0(:,2);
[Nele0,nod0]=size(connectivity0);
Nnodes0=max(connectivity0(:));

%% Error estimate of solution on a given grid

%% Implicit error estimate
%Subdivide all elements and perform a diagnostic soluion on that grid

[coordinates1,connectivity1]=FE2dRefineMesh(coordinates0,connectivity0);



x1=coordinates1(:,1) ; y1=coordinates1(:,2);
Nnodes1=max(connectivity1(:)); [Nele1,nod1]=size(connectivity1);
fprintf(' Error estimate meshrefinement %-2.2i \t # Elements: %-i \t \t # Nodes=%-i \n',iRefinement,Nele1,Nnodes1)
% Calculate solution on grid1 (which is twice as fine as grid0) and used the solution to identify which elements of Grid0 need to be subdivided
[DTxy1,TRIxy1,u1,v1,h1,s1,b1,S1,B1,AGlen1,C1,rho1]=DiagnosticSolutionGrid1toGrid2(Experiment,coordinates1,connectivity1,...
    Boundary0,MeshBoundaryCoordinates,u0,v0,h0,s0,b0,AGlen0,C0,n,m,rho,rhow,g,Itime,nip,DTxy0,...
    CtrlVar);


I=1; % this is a cell array, and possible a number of different refinementcriterie are defnie, however here
     % currently only the first one is considerd
fprintf(' error estimate based on : %s \n ',CtrlVar.RefineCriteria{I})

switch lower(CtrlVar.RefineCriteria{I})
    
    case 'effective strain rates';
        
        nip=1;
        save TestSave
        [~,xint1,yint1,exx1,eyy1,exy1,~,e1]=calcStrainRatesEtaInt(u1,v1,coordinates1,connectivity1,nip,AGlen1,n,CtrlVar);
        [~,xint0,yint0,exx0,eyy0,exy0,~,e0]=calcStrainRatesEtaInt(u0,v0,coordinates0,connectivity0,nip,AGlen0,n,CtrlVar);
         
        DTint1 = DelaunayTri(xint1,yint1);
        exx1=Grid1toGrid2(DTint1,exx1,xint0,yint0);
        eyy1=Grid1toGrid2(DTint1,eyy1,xint0,yint0);
        exy1=Grid1toGrid2(DTint1,exy1,xint0,yint0);
        
        eError0=real(sqrt((exx0-exx1).^2+(eyy0-eyy1).^2+(exx0-exx1).*(eyy0-eyy1)+(exy0-exy1).^2)); % ele values

        EleError=eError0;
        NormalizedEleError=EleError/max(e0);
        
        %figure ; trisurf(TRIxy0,x0,y0,e0) ;  title(sprintf(' effective strain rate on nodes for grid0 I=%-i ',iRefinement))
        %figure ; trisurf(TRIxy1,x1,y1,e1) ;  title(sprintf(' effective strain rate on nodes for grid1 I=%-i ',iRefinement))
        
        %figure
        %patch('Faces',connectivity0,'Vertices',coordinates0,'FaceVertexCData',eError0,'CDataMapping','scaled','FaceColor','flat');
        %title(sprintf(' eError0 on grid0 I=%-i ',iRefinement)) ; colorbar ; axis equal
       
        
        
    case 'speed'
        
        % compare solhations on original grid and the finer grid, and select which elements to subdivide
        speed1=sqrt(u1.*u1+v1.*v1);
        speed0=sqrt(u0.*u0+v0.*v0);
        u1Grid0=Grid1toGrid2(DTxy1,u1,x0,y0); v1Grid0=Grid1toGrid2(DTxy1,v1,x0,y0);
        
        
        figure ; trisurf(TRIxy0,x0,y0,speed0) ; title(sprintf(' speed0 on grid0 I=%-i Nele=%-i',iRefinement,Nele0))
        figure ; trisurf(TRIxy1,x1,y1,speed1) ; title(sprintf(' speed1 on grid1 I=%-i Nele=%-i',iRefinement,Nele1))
        
        SpeedError0=sqrt((u0-u1Grid0).*(u0-u1Grid0)+(v0-v1Grid0).*(v0-v1Grid0));
        
        figure ; trisurf(TRIxy0,x0,y0,SpeedError0) ;  title(sprintf(' SpeedError on grid0 I=%-i ',iRefinement))
        
        SpeedError0EleAverage=real(mean(reshape(SpeedError0(connectivity0),Nele0,nod0),2));
        M0 = Ele2Nodes(connectivity0,Nnodes0);
        %SpeedError0=M0*SpeedError0EleAverage;
        SpeedError0=SpeedError0EleAverage;
        NormalizedSpeedError0=SpeedError0/max(speed0);
        
        figure ; trisurf(TRIxy0,x0,y0,M0*NormalizedSpeedError0) ;  title(sprintf(' NormalizedSpeedError0 averaged on grid0 I=%-i ',iRefinement))
        
        %%
        EleError=SpeedError0;
        NormalizedEleError=EleError/max(speed0);
        
    otherwise
        error(' what case? ')
end

rmsEleError=norm(EleError)/Nele0;


%figure
%hist(EleError,20)
%title(sprintf(' eError0 on grid0 I=%-i ',iRefinement)) ; 

return

end