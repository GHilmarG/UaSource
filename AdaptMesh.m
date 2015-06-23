
function [CtrlVar,MUAnew,BCsNew,MeshBoundaryCoordinates,GF,GLdescriptors,...
    s,b,h,S,B,ub,vb,ud,vd,ubvbLambda,udvdLambda,rho,rhow,g,AGlen,n,C,m,ab,as,dhdt,dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1]=...
    AdaptMesh(CtrlVar,Experiment,MeshBoundaryCoordinates,MUAold,BCsOld,time,Itime,...
    GF,GLdescriptors,alpha,...
    s,b,h,S,B,ub,vb,ud,vd,Ruv,Lubvb,ubvbLambda,udvdLambda,rho,rhow,g,AGlen,n,C,m,ab,as,dhdt,dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1)

narginchk(42,42)
persistent MUA_Background

%%
%save TestSave
%error('fsda')




MUAnew=MUAold;
hOld=h;
BCsNew=BCsOld;
MeshAdvanceRetreat=0;
CtrlVar.MeshChanged=0;

% I want to control directly when FE-mesh is plotted, and avoid the potential FE-mesh plotting in genmesh2d
PlotMeshOnInput=CtrlVar.PlotMesh; CtrlVar.PlotMesh=0;


if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=2
    figure
    PlotFEmesh(MUAold.coordinates,MUAold.connectivity,CtrlVar);
    title(sprintf('Initial mesh on input to AdaptMesh : #Ele=%-i, #Nodes=%-i, #nod=%-i',MUAold.Nele,MUAold.Nnodes,MUAold.nod))
end



if CtrlVar.TimeGeometries.Flag
    
    %%  global remeshing due to imposed changes in geometry, e.g. calving event
    CtrlVar.MeshChanged=0;
    N=numel(CtrlVar.TimeGeometries.StartTime)  ;  % Total number of possible geometry changes with time
    M=CtrlVar.TimeGeometries.Done              ;  % The number of `processed'geometry changes
    if M < N    % still some more changes in geometry to be processed
        
        if time >= CtrlVar.TimeGeometries.StartTime(M+1)
            
            if CtrlVar.InfoLevelAdaptiveMeshing>=1
                fprintf(CtrlVar.fidlog,' Geometry changed by user, global remeshing at time %-g \n ',time);
            end
            
            load(CtrlVar.TimeGeometries.MeshBoundaryInputFile{M+1},'MeshBoundaryCoordinates')
            CtrlVar.GmeshFile=CtrlVar.TimeGeometries.GmeshFile{M+1};
            CtrlVar.GmeshMeshingMode=CtrlVar.TimeGeometries.GmeshMeshingMode{M+1};
            CtrlVar.TimeGeometries.Done=M+1;
            if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=10
                figure ; PlotFEmesh(MUAold.coordinates,MUAold.connectivity,CtrlVar);
                title(sprintf('Before time-geometry %-i. # Ele=%-i, #Nodes=%-i, #nod=%-i',...
                    CtrlVar.TimeGeometries.Done,MUAold.Nele,MUAold.Nnodes,MUAold.nod))
            end
            
            
            MUAnew=genmesh2d(CtrlVar,MeshBoundaryCoordinates);
            
            if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=10
                figure ; PlotFEmesh(MUAnew.coordinates,MUAnew.connectivity,CtrlVar);
                title(sprintf('After time-geometry %-i. # Ele=%-i, #Nodes=%-i, #nod=%-i',...
                    CtrlVar.TimeGeometries.Done,MUAnew.Nele,MUAnew.Nnodes,MUAnew.nod))
            end
            
            CtrlVar.MeshChanged=1;
        end
        
    end
    
elseif CtrlVar.FEmeshAdvanceRetreat && ( ReminderFraction(time,CtrlVar.FEmeshAdvanceRetreatDT)<1e-5 || CtrlVar.FEmeshAdvanceRetreatDT==0)
    
    MeshAdvanceRetreat=1;
    %%  FE mesh advance or retreat
    if CtrlVar.InfoLevelAdaptiveMeshing>=1
        fprintf(CtrlVar.fidlog,' FE mesh advance or retreat, element activated/deactivated at time %-g \n ',time);
    end
    
    if isempty(MUA_Background)
        
        if    exist(fullfile(cd,CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName), 'file')  == 2 ...
                || exist(fullfile(cd,[CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName,'.mat']), 'file')  == 2
            try
                
                if CtrlVar.InfoLevelAdaptiveMeshing>=1
                    fprintf('Reading ''MUA_Background'' from: %s ',CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName)
                end
                
                %load(CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName,'coordinates','connectivity')
                load(CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName,'MUA_Background')
                
                if exist('MUA_Background','var')==0
                    fprintf(' The variable MUA_Background not found in %s.! \n',CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName)
                    error('AdaptMesh:MUA_BackgroundNotFound','Where is MUA_Background?')
                end
                
                MUA_Background=UpdateMUA(CtrlVar,MUA_Background);
                if CtrlVar.InfoLevelAdaptiveMeshing>=1
                    fprintf('done \n ')
                end
               
            catch
                fprintf('File %s not containing ''MUA_Background'' \n',CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName)
                error(' exiting ' )
            end
        else
            fprintf('File with background meshfile %s could not be found \n',CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName)
            error(' exiting ' )
        end

    end
    
    % map thickness onto backgroundmeshfile
    
    % I should replace this with a `mapping' from the existing mesh onto the background mesh,
    % rather than use this general interpolation routine.
    
    %%
    %hBackground=Grid1toGrid2(DTxy,h,coordinatesBackground(:,1),coordinatesBackground(:,2),CtrlVar,CtrlVar.ThickMin);
    hBackground=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,MUA_Background.coordinates(:,1),MUA_Background.coordinates(:,2),CtrlVar.ThickMin,h);
    iDeactivatedElements=FindElementsToDeactivate(CtrlVar,MUA_Background,hBackground);
    
    fprintf(CtrlVar.fidlog,'%i elements of background mesh deactivated. ',numel(find(iDeactivatedElements)));
    
    
    if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=10
        
        figure(10000)
        subplot(2,1,1,'replace')
        hold off
        PlotNodalBasedQuantities(MUA_Background.connectivity,MUA_Background.coordinates,hBackground,CtrlVar);
        xlabel('x (km)') ; ylabel('y (km)') ;
        colorbar ; title(colorbar,'(m)') ; % caxis([0 100])
        hold on
        
        title(sprintf('hBackground at t=%#5.1f ',time))
        %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
        axis equal tight
        
        ind=hBackground<=CtrlVar.ThickMinDeactivateElements;
        plot(MUA_Background.coordinates(ind,1)/CtrlVar.PlotXYscale,MUA_Background.coordinates(ind,2)/CtrlVar.PlotXYscale,'.r')
        hold off
        
        
        subplot(2,1,2,'replace')
        hold off
        PlotFEmesh(MUA_Background.coordinates,MUA_Background.connectivity,CtrlVar)
        hold on
        ind=hBackground<=CtrlVar.ThickMinDeactivateElements;
        plot(MUA_Background.coordinates(ind,1)/CtrlVar.PlotXYscale,MUA_Background.coordinates(ind,2)/CtrlVar.PlotXYscale,'.r')
        
        xEleBackground=Nodes2EleMean(MUA_Background.connectivity,MUA_Background.coordinates(:,1));
        yEleBackground=Nodes2EleMean(MUA_Background.connectivity,MUA_Background.coordinates(:,2));
        
        plot(xEleBackground(iDeactivatedElements)/CtrlVar.PlotXYscale,yEleBackground(iDeactivatedElements)/CtrlVar.PlotXYscale,'og')
        title(sprintf('time=%-g \t #Ele=%-i, #Nodes=%-i, #nod=%-i, # EleBG=%-i, #NodesBG=%-i, #nodBG=%-i',...
            CtrlVar.time,MUAold.Nele,MUAold.Nnodes,MUAold.nod,MUA_Background.Nele,MUA_Background.nod,MUA_Background.Nnodes))
        axis equal tight
    end
    
    % eliminate elements
    if numel(find(iDeactivatedElements))>0
        [coordinates,connectivity]=DeactivateElements(CtrlVar,iDeactivatedElements,MUA_Background.coordinates,MUA_Background.connectivity);
        clear MUA
        MUAnew=CreateMUA(CtrlVar,connectivity,coordinates);
        clear connectivity coordinates
        MUAnew=UpdateMUA(CtrlVar,MUAnew);
        
        CtrlVar.MeshChanged=1;
        fprintf(CtrlVar.fidlog,'Change in the number of elements : %+i \n ',MUAnew.Nele-MUAold.Nele);
    end
    
end


%%

iAdapt=CtrlVar.AdaptMesh  ...
    && (CtrlVar.AdaptMeshInitial || (mod(Itime,CtrlVar.AdaptMeshInterval)==0 && CtrlVar.AdaptMeshInterval>0 && Itime>1)) ...
    && ~MeshAdvanceRetreat;



% if mesh has changed, for example in the MeshAdvanceRetreat step, then I continue
if ~iAdapt && ~CtrlVar.MeshMorphing && ~CtrlVar.MeshChanged
    CtrlVar.PlotMesh=PlotMeshOnInput;
    return
end


if CtrlVar.MeshMorphing || CtrlVar.TimeGeometries.Flag  || MeshAdvanceRetreat % no iterations needed for mesh morphing and TimeGeometries
    Iterations=1;
else
    Iterations=CtrlVar.AdaptMeshIterations;
end

if CtrlVar.InfoLevelAdaptiveMeshing>=1
    fprintf(CtrlVar.fidlog,'Before remeshing: '); PrintInfoAboutElementsSizes(CtrlVar,MUAold)
end

for JJ=1:Iterations
    if iAdapt
        
        CtrlVar.AdaptMeshInitial=0;
        if CtrlVar.InfoLevelAdaptiveMeshing>=1
            fprintf(CtrlVar.fidlog,' =====  Remeshing at start of time step %-i. Iteration #%-i out of %-i \n ',Itime,JJ,Iterations);
        end
        
        switch lower(CtrlVar.MeshRefinementMethod)
            
            case 'implicit'
                
                error(' implicit error estimate not fully implemented \n')
                [coordinates,connectivity]=DesiredEleSizesBasedOnImplicitErrorEstimate(Experiment,MeshBoundaryCoordinates,Boundary,...
                    s,b,S,B,h,ub,vb,coordinates,connectivity,nip,AGlen,C,Luv,Luvrhs,ubvbLambda,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
                
            case {'explicit:global','explicit:local'}
                
                if Itime==1
                    
                    % do I need to calculate velocities for error estimates?
                    CalcVel=0;
                    for J=1:length(CtrlVar.RefineCriteria)
                        if strcmp(CtrlVar.RefineCriteria{J},'effective strain rates') ; CalcVel=1 ; end
                        if strcmp(CtrlVar.RefineCriteria{J},'residuals') ; CalcVel=1 ; end
                        if CalcVel
                            break
                        end
                    end
                    
                    if CalcVel
                        MUAold=UpdateMUA(CtrlVar,MUAold);
                         [ub,vb,ud,vd,ubvbLambda,udvdLambda,Kuv,Ruv,RunInfo,Lubvb]= uv(CtrlVar,MUAold,BCsOld,s,b,h,S,B,ub,vb,ud,vd,ubvbLambda,udvdLambda,AGlen,C,n,m,alpha,rho,rhow,g,GF);
                    end
 
                end
                
                %save TestSave ; error('afds')
                [MUAnew,xGLmesh,yGLmesh,CtrlVar]=...
                    RemeshingBasedOnExplicitErrorEstimate(MeshBoundaryCoordinates,...
                    S,B,h,s,b,ub,vb,dhdt,MUAold,AGlen,C,n,rho,rhow,CtrlVar,GF,Ruv,Lubvb,ubvbLambda);
                
                CtrlVar.MeshChanged=1;
                
            otherwise
                error(' unknown case  ')
        end
        
        % I might be doing some GL Mesh Morphing at a later stage
        % if so, then I need to recalculate the GLdescriptors
        % This allows for global remeshing to be done once in a while
        % with frequent GL Mesh Morphing in between.
        if CtrlVar.GLmeshing && CtrlVar.MeshMorphing
            
            [dtGL,GLdescriptors,xGLbackground,yGLbackground]=...
                CreateBackgroundGLmesh(MUAnew.coordinates,MeshBoundaryCoordinates,xGLmesh,yGLmesh,CtrlVar);
            
             if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=10
                figure ; PlotFEmesh(MUAnew.coordinates,MUAnew.connectivity,CtrlVar);
                hold on ; plot(xGLbackground,yGLbackground,'-go','LineWidth',2);
                hold on ; plot(xGLmesh,yGLmesh,'-rx','LineWidth',2);
                title('New GL background mesh')
            end
        end
        
    elseif CtrlVar.MeshMorphing && ~isempty(dtGL)
        
        if CtrlVar.InfoLevelAdaptiveMeshing>=1
            fprintf(CtrlVar.fidlog,'GLmorphing \n');
        end
        
        [coordinates,xGLmorphing,yGLmorphing]=GLmorphing(CtrlVar,Experiment,MUAold.coordinates,MUAold.connectivity,GF,MeshBoundaryCoordinates,GLdescriptors,dtGL);
        connectivity=MUAold.connectivity;
        MUAnew=CreateMUA(CtrlVar,connectivity,coordinates);
        
        CtrlVar.MeshChanged=1;
        notOK=1; maxit=30; tol=0.5;
        if CtrlVar.InfoLevelAdaptiveMeshing>=1 ; fprintf(CtrlVar.fidlog,'After GL morphing '); end
        while notOK
            q = MeshQuality(coordinates,connectivity);
            Q=0.1;
            iq=q<Q; notOK=any(iq) ;
            if notOK
                fprintf(CtrlVar.fidlog,'%i elements with quality factor < %f \n',sum(iq),Q);
                fprintf(CtrlVar.fidlog,' Laplace smoothing of mesh \n');
                [MUAnew.coordinates,MUAnew.connectivity]=FEmeshSmoothing(MUAnew.coordinates,MUAnew.connectivity,maxit,tol);
            else
                fprintf(CtrlVar.fidlog,'all elements with quality factor > %f \n',Q);
            end
            tol=tol/2;
            
        end
        
        if CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=10
            %%
            figure ; PlotFEmesh(MUAnew.coordinates,MUAnew.connectivity,CtrlVar);
            hold on ; plot(xGLmorphing,yGLmorphing,'-gx','LineWidth',2);
            title('Mesh after GLmorphing')
            %%
        end
    end
    
    
    %%
    % no further modification of coordinates or connectivity
    %
    CtrlVar.MeshChanged=1;
    
    
    if CtrlVar.InfoLevelAdaptiveMeshing>=1;
        fprintf(CtrlVar.fidlog,'After remeshing: ') ; PrintInfoAboutElementsSizes(CtrlVar,MUAnew)
    end
    
    if MUAnew.Nele==0
        fprintf('No elements left in mesh! \n ')
        CtrlVar.PlotMesh=PlotMeshOnInput;
        return
    end
     
    OutsideValues=[0 ; 0 ; 0; 0 ; 0; 0 ; 0 ; 0 ; 0 ; 0; 0 ; 0 ; 0 ; 0 ];
    [s,b,h,S,B,rho,AGlen,n,C,m,GF,ub,vb,ud,vd,dhdt,dubdt,dvbdt,duddt,dvddt,dhdtm1,dubdtm1,dvbdtm1,duddtm1,dvddtm1]=...
        MapQuantitiesToNewFEmesh(CtrlVar,MUAnew,MUAold,hOld,time,OutsideValues,...
        ub,vb,ud,vd,dhdt,dubdt,dvbdt,duddt,dvddt,dhdtm1,dubdtm1,dvbdtm1,duddtm1,dvddtm1);
    
    BCsNew=BoundaryConditions;
    BCsNew=GetBoundaryConditions(Experiment,CtrlVar,MUAnew,BCsNew,time,s,b,h,S,B,ub,vb,ud,vd,GF);
        
    [as,ab]=DefineMassBalance(Experiment,CtrlVar,MUAnew,time,s,b,h,S,B,rho,rhow,GF);
    
    if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=1
        if CtrlVar.PlotBCs
            figure ; PlotBoundaryConditions(CtrlVar,MUAnew,BCsNew);
        else
            figure
            PlotFEmesh(MUAnew.coordinates,MUAnew.connectivity,CtrlVar); hold on
            title(sprintf('AdaptMesh Iteration %i: #Ele=%-i, #Nodes=%-i, #nod=%-i',JJ,MUAnew.Nele,MUAnew.Nnodes,MUAnew.nod))
        end
        
        if CtrlVar.GLmeshing
            hold on ; plot(xGLmesh,yGLmesh,'-rx','LineWidth',2);
        end
        
    end
      
    if any(isnan(C)) ; error( ' C nan ') ; end
    if any(isnan(AGlen)) ; error( ' AGlen nan ') ; end
    if any(isnan(S)) ; error( ' S nan ') ; end
    if any(isnan(h)) ; error( ' h nan ') ; end
    if any(isnan(ub)) ; error( ' u nan ') ; end
    if any(isnan(vb)) ; error( ' v nan ') ; end
    if any(isnan(rho)) ; error( ' rho nan ') ; end
    %%
    
  
    
    if JJ<Iterations
        
        % do I need to calculate velocities for error estimates?
        CalcVel=0;
        for J=1:length(CtrlVar.RefineCriteria)
            if strcmp(CtrlVar.RefineCriteria{J},'effective strain rates') ; CalcVel=1 ; end
        end
        
        if CalcVel
            
            ub=ub*0 ; vb=vb*0 ; ubvbLambda=ubvbLambda*0; % experience has shown that it is almost always best here to reset estimates of (u,v) to zero
            ud=ud*0 ; vd=vd*0;
            MUAnew=UpdateMUA(CtrlVar,MUAnew);
            [ub,vb,ud,vd,ubvbLambda,udvdLambda,Kuv,Ruv,RunInfo]= uv(CtrlVar,MUAnew,BCsNew,s,b,h,S,B,ub,vb,ud,vd,ubvbLambda,udvdLambda,AGlen,C,n,m,alpha,rho,rhow,g,GF);
        end
        
        
        
    end
    MUAnew=UpdateMUA(CtrlVar,MUAnew);
    MUAold=MUAnew;  hOld=h; BCsOld=BCsNew;
end

%MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);
MUAnew=UpdateMUA(CtrlVar,MUAnew);

CtrlVar.InitialDiagnosticStep=1;  % make sure that in next uvh step I start with an initial uv step

if ~isempty(CtrlVar.SaveAdaptMeshFileName)
    MUA=MUAnew;
    save(CtrlVar.SaveAdaptMeshFileName,'MUA')
    fprintf(CtrlVar.fidlog,' Adapted FE mesh was saved in %s .\n',CtrlVar.SaveAdaptMeshFileName);
end

CtrlVar.PlotMesh=PlotMeshOnInput;

end

