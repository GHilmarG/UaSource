function [MUA,xGLmesh,yGLmesh,CtrlVar]=...
    RemeshingBasedOnExplicitErrorEstimate(MeshBoundaryCoordinates,S0,B0,h0,s0,b0,ub0,vb0,ud0,vd0,dhdt0,MUA,AGlen0,C0,n,rho0,rhow,CtrlVar,GF0,Ruv,Lubvb,ubvbLambda)


%save TestSave ; error('dfsa')
%

% Global remeshing, mesh morhing, or local mesh refinement based on explicit error estimate.
%
% On output MUA holds the new mesh. On input MUA is an existing mesh
% xGLmesh and yGLmesh is only relevant for GL morphing
%

xGLmesh=[] ; yGLmesh=[];


S=S0 ; B=B0;h=h0; rho=rho0;
%u=u0;v=v0 ; dhdt=dhdt0; AGlen=AGlen0;


hf=(S-B)*rhow./rho ;

%%    Step 1 : Define desired size of elements based on some criteria
% x, y are x,y coordinates of nodes
[x0,y0,EleSizeDesired,EleSizeCurrent,ElementsToBeRefined,NodalErrorIndicators]=DesiredEleSizes(CtrlVar,MUA,s0,b0,S0,B0,rho0,rhow,ub0,vb0,ud0,vd0,dhdt0,h0,hf,AGlen0,n,GF0,Ruv,Lubvb,ubvbLambda);

if strcmp(CtrlVar.MeshGenerator,'gmsh')
    if norm([x0-MUA.coordinates(:,1);y0-MUA.coordinates(:,2)]) > 100*eps
        error('RemeshingBasedOnExplicitErrorEstimate:gmsh','When using gmsh desired ele sizes must be defined at nodes')
    end
end


TRIxy0=TriFE(MUA.connectivity);  % this is the triangulation of the input FEmesh over which
% the error estimation is performed

if any(isnan(EleSizeDesired)) ; save TestSave ; error('fdsa') ; end

%% mesh refinement or global remeshing


switch lower(CtrlVar.MeshRefinementMethod)
    
    case 'explicit:local'
        

        
        MUA=LocalMeshRefinement(CtrlVar,MUA,ElementsToBeRefined);
        
%         % refine and smoothmesh only works for 3-nod elements
%         [MUA.coordinates,MUA.connectivity]=ChangeElementType(MUA.coordinates,MUA.connectivity,3);
%         [MUA.coordinates,MUA.connectivity] = refine(MUA.coordinates,MUA.connectivity,ElementsToBeRefined);
%         [MUA.coordinates] = GHGsmoothmesh(MUA.coordinates,MUA.connectivity,CtrlVar.LocalAdaptMeshSmoothingIterations,[]);
%          MUA.connectivity=FlipElements(MUA.connectivity); 
%           
%         % create MUA (this takes care of any eventual change in element type as well)
%         MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);

           
    case 'explicit:global'
        %% Global remeshing
        switch lower(CtrlVar.MeshGenerator)
            case 'mesh2d'
                CtrlVar.MeshSize=zeros(length(x0),3);
                CtrlVar.MeshSize(:,1)=x0 ; CtrlVar.MeshSize(:,2)=y0; CtrlVar.MeshSize(:,3)=EleSizeDesired;
            case 'gmsh'
                
                GmshBackgroundScalarField.xy=[x0(:) y0];
                GmshBackgroundScalarField.EleSize=EleSizeDesired(:) ;
                GmshBackgroundScalarField.TRI=TRIxy0;
                
                
                %                 if   CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1
                %                     figure ; PlotFEmesh(GmshBackgroundScalarField.xy,GmshBackgroundScalarField.TRI,CtrlVar);
                %                 end
                
            otherwise
                error('Mesh generator not correctly defined. Define variable CtrlVar.MeshGenerator {mesh2d|gmsh} ')
        end
        
        if CtrlVar.GLmeshing==1
            GF = GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);
            [MeshBoundaryCooWithGLcoo,edge,face,xGLmesh,yGLmesh]=glLineEdgesFaces(GF,MUA.coordinates,MUA.connectivity,MeshBoundaryCoordinates,CtrlVar);
            MUA=genmesh2d(CtrlVar,MeshBoundaryCooWithGLcoo,edge,face);
        else
            MUA=genmesh2d(CtrlVar,MeshBoundaryCoordinates,[],[],GmshBackgroundScalarField);
        end
        
        %figure ; PlotFEmesh(coordinates,connectivity,CtrlVar);
        
        
        
        % if the ratio of actual number of elements to desired number of elements
        % is either larger than `EleFactorU' or smaller than 'EleFactorDown' then MeshSizeMin is scaled
        % up or down by the factor 'EleFactor' and a new remeshing is done.
        % The maximum size of elements does not change and is always MeshSizeMax. Only MeshSizeMin is changed
        % together with the error chriteria, MeshSizeMax and the number of elements are the main variables,
        % affecting the mesh. The new value of MeshSizeMin is written out and in principle one should use that
        % value of MeshSizeMin as an input value in future runs of the same model
        
        %EleFactorUp=1.3 ; EleFactorDown=0.1;
        EleFactorUp=CtrlVar.MaxNumberOfElementsUpperLimitFactor;
        EleFactorDown=CtrlVar.MaxNumberOfElementsLowerLimitFactor;
        
        It=0;
        
        while (MUA.Nele>EleFactorUp*CtrlVar.MaxNumberOfElements ||  MUA.Nele<EleFactorDown*CtrlVar.MaxNumberOfElements ) && It<4
            
            ScalingFactor=sqrt(MUA.Nele/CtrlVar.MaxNumberOfElements);
            
            %EleSize=EleSize*ScalingFactor;
            % make sure that no element is greater than max element size
            %CtrlVar.MeshSize(EleSize>CtrlVar.MeshSizeMax,3)=CtrlVar.MeshSizeMax;
            
            % rescale elements sizes so that maximum size is still MeshSizeMax but min ele size is scaled either up or down
            maxE=max(EleSizeDesired);
            minE=min(EleSizeDesired);
            
            NewMinE=min([minE*ScalingFactor,0.9*CtrlVar.MeshSizeMax]);
            if maxE~=minE
                EleSizeDesired=NewMinE+(CtrlVar.MeshSizeMax-NewMinE)*(EleSizeDesired-minE)/(maxE-minE);
            else
                EleSizeDesired=NewMinE+EleSizeDesired*0;
            end
            
            It=It+1;
            
            if NewMinE ~=CtrlVar.MeshSizeMin
                warning('RmeshingBasedOnExplicitErrorEstimate:NewMeshSizeMin','CtrlVar.MeshSizeMin set to %f',CtrlVar.MeshSizeMin)
            end
            
            CtrlVar.MeshSizeMin=NewMinE;
            
            
            switch lower(CtrlVar.MeshGenerator)
                case 'mesh2d'
                    CtrlVar.MeshSize=zeros(length(x0),3);
                    CtrlVar.MeshSize(:,1)=x0 ; CtrlVar.MeshSize(:,2)=y0; CtrlVar.MeshSize(:,3)=EleSizeDesired;
                case 'gmsh'
                    GmshBackgroundScalarField.xy=[x0(:) y0(:)] ;
                    GmshBackgroundScalarField.EleSize=EleSizeDesired(:) ;
                    GmshBackgroundScalarField.TRI=TRIxy0 ;
                    
                    if any(isnan(EleSizeDesired)) ; error('fdsa') ; end
                otherwise
                    error('Mesh generator not correctly defined. Define variable CtrlVar.MeshGenerator {mesh2d|gmsh} ')
            end
            
            if CtrlVar.InfoLevelAdaptiveMeshing>=1;
                if MUA.Nele>EleFactorUp*CtrlVar.MaxNumberOfElements
                    fprintf(CtrlVar.fidlog,'Nele=%-i > %-i , hence desired meshsize is scaled up. New MeshSizeMin is %-g \n',MUA.Nele,CtrlVar.MaxNumberOfElements,CtrlVar.MeshSizeMin);
                else
                    fprintf(CtrlVar.fidlog,'Nele=%-i < %-i , hence desired meshsize is scaled down. New MeshSizeMin is %-g \n',MUA.Nele,CtrlVar.MaxNumberOfElements,CtrlVar.MeshSizeMin);
                end
            end
            
            if CtrlVar.GLmeshing==1
                MUA=genmesh2d(CtrlVar,MeshBoundaryCooWithGLcoo,edge,face);
            else
                MUA=genmesh2d(CtrlVar,MeshBoundaryCoordinates,[],[],GmshBackgroundScalarField);
            end
            
            if CtrlVar.InfoLevelAdaptiveMeshing>=1;
                fprintf(CtrlVar.fidlog,'new Nele after rescaling is %-i and CtrlVar.MaxNumberOfElements is %-i  \n',MUA.Nele,CtrlVar.MaxNumberOfElements);
            end
        end
        
    otherwise
        fprintf('Incorrect value for CtrlVar.MeshRefinementMethod (%s) \n',CtrlVar.MeshRefinementMethod)
        error('RemeshingBasedOnExplicitErrorEstimate:CaseNotFound','Case not found.')
end



end

