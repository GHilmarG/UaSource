function [coordinates,connectivity,xGLmesh,yGLmesh,CtrlVar]=...
    DesiredEleSizesBasedOnExplicitErrorEstimate(Experiment,MeshBoundaryCoordinates,...
    S0,B0,h0,s0,b0,u0,v0,wSurf0,dhdt0,coordinates,connectivity,nip,AGlen0,C0,n,rho0,rhow,DTxy0,TRIxy0,CtrlVar)

xGLmesh=[] ; yGLmesh=[];

Itime=[] ; time=[];



x=coordinates(:,1); y=coordinates(:,2);
%if I==1
    S=S0 ; B=B0;h=h0;u=u0;v=v0;wSurf=wSurf0 ; dhdt=dhdt0; AGlen=AGlen0; rho=rho0; DTxy=DTxy0; TRIxy=TRIxy0;

hf=(S-B)*rhow./rho ;

%    Step 1 : Define desired size of elements based on some criteria
[xEle,yEle,EleSize]=DesiredEleSizes(u,v,wSurf,dhdt,h,hf,coordinates,connectivity,nip,AGlen,n,DTxy,TRIxy,CtrlVar);

if any(isnan(EleSize)) ; error('fdsa') ; end

if   CtrlVar.doplots==1 && CtrlVar.doRemeshPlots==1; figure(123) ; trisurf(TRIxy,xEle,yEle,EleSize) ;  title(' desired ele sizes ') ; end


if CtrlVar.MeshRefinement==1
    
    assert(size(connectivity,2)==3,' Local mesh refinement currently only implemented for 3Nod elements')
    RefineNodes=double(EleSize < (min(EleSize)+0.5*( max(EleSize)-min(EleSize))));
    trisurf(TRIxy,xEle,yEle,RefineNodes) ;  title(' RefineNodes')
    
    RefineElements=sum(RefineNodes(connectivity),2) > 2 ; % if more than two nodes
    
    [coordinates,connectivity] = refine(coordinates,connectivity,RefineElements);
    [coordinates,connectivity] = smoothmesh(coordinates,connectivity);
    if   CtrlVar.doplots==1 && CtrlVar.doRemeshPlots==1
        figure
        patch('faces',connectivity,'vertices',coordinates,'facecolor','none','edgecolor','b');
    end
    
else
    %% Global remeshing
    switch lower(CtrlVar.MeshGenerator)
        case 'mesh2d'
            CtrlVar.MeshSize=zeros(length(xEle),3);
            CtrlVar.MeshSize(:,1)=xEle ; CtrlVar.MeshSize(:,2)=yEle; CtrlVar.MeshSize(:,3)=EleSize;
        case 'gmesh'
            
            GmeshBackgroundScalarField.xy=[xEle(:) yEle];
            GmeshBackgroundScalarField.EleSize=EleSize(:) ;
            GmeshBackgroundScalarField.TRI=TRIxy ;
            
            if   CtrlVar.doplots==1 && CtrlVar.doRemeshPlots==1
                figure ; PlotFEmesh(GmeshBackgroundScalarField.xy,GmeshBackgroundScalarField.TRI,CtrlVar);
            end
            
        otherwise
            error('Mesh generator not correctly defined. Define variable CtrlVar.MeshGenerator {mesh2d|gmesh} ')
    end
    
    if CtrlVar.GLmeshing==1
        GF = GL2d(B,S,h,rhow,rho,connectivity,CtrlVar);
        [MeshBoundaryCooWithGLcoo,edge,face,xGLmesh,yGLmesh]=glLineEdgesFaces(GF,coordinates,connectivity,MeshBoundaryCoordinates,CtrlVar);
        [coordinates,connectivity]=genmesh2d(Experiment,MeshBoundaryCooWithGLcoo,CtrlVar,edge,face);
    else
        [coordinates,connectivity]=genmesh2d(Experiment,MeshBoundaryCoordinates,CtrlVar,[],[],GmeshBackgroundScalarField);
    end
    
    %figure ; PlotFEmesh(coordinates,connectivity,CtrlVar);
    [Nele,~]=size(connectivity);
    
    
    % if the ratio of actual number of elements to desired nbumber of elements
    % is either larger than `EleFactorU' or smaller than 'EleFactorDown' then MeshSizeMin is scaled
    % up or down by the factor 'EleFactor' and a new remeshing is done.
    % The maximum size of elements does not change and is always MeshSizeMax. Only MeshSizeMin is changed
    % together with the error chriteria, MeshSizeMax and the number of elements are the main variables,
    % affecting the mesh. The new value of MeshSizeMin is written out and in principle one should use that
    % value of MeshSizeMin as an input value in future runs of the same model
    
    EleFactorUp=1.3 ; EleFactorDown=0.1; It=0;
    while (Nele>EleFactorUp*CtrlVar.MaxNumberOfElements ||  Nele<EleFactorDown*CtrlVar.MaxNumberOfElements ) && It<4
        
        ScalingFactor=sqrt(Nele/CtrlVar.MaxNumberOfElements);
        
        %EleSize=EleSize*ScalingFactor;
        % make sure that no element is greater than max element size
        %CtrlVar.MeshSize(EleSize>CtrlVar.MeshSizeMax,3)=CtrlVar.MeshSizeMax;
        
        % rescale elements sizes so that maximum size is still MeshSizeMax but min ele size is scaled either up or down
        maxE=max(EleSize);
        minE=min(EleSize);
        fprintf(' maxE=%-g \t minE=%-g \n ',maxE,minE)
        NewMinE=min([minE*ScalingFactor,0.9*CtrlVar.MeshSizeMax]);
        EleSize=NewMinE+(CtrlVar.MeshSizeMax-NewMinE)*(EleSize-minE)/(maxE-minE);
        
        It=It+1;
        CtrlVar.MeshSizeMin=NewMinE;
        
        
        switch lower(CtrlVar.MeshGenerator)
            case 'mesh2d'
                CtrlVar.MeshSize=zeros(length(xEle),3);
                CtrlVar.MeshSize(:,1)=xEle ; CtrlVar.MeshSize(:,2)=yEle; CtrlVar.MeshSize(:,3)=EleSize;
            case 'gmesh'
                GmeshBackgroundScalarField.xy=[xEle(:) yEle(:)] ;
                GmeshBackgroundScalarField.EleSize=EleSize(:) ;
                GmeshBackgroundScalarField.TRI=TRIxy ;
                if any(isnan(EleSize)) ; error('fdsa') ; end
            otherwise
                error('Mesh generator not correctly defined. Define variable CtrlVar.MeshGenerator {mesh2d|gmesh} ')
        end
        
        if Nele>EleFactorUp*CtrlVar.MaxNumberOfElements
            fprintf(CtrlVar.fidlog,'Nele=%-i > %-i , hence desired meshsize is scaled up. New MeshSizeMin is %-g \n',Nele,CtrlVar.MaxNumberOfElements,CtrlVar.MeshSizeMin);
        else
            fprintf(CtrlVar.fidlog,'Nele=%-i < %-i , hence desired meshsize is scaled down. New MeshSizeMin is %-g \n',Nele,CtrlVar.MaxNumberOfElements,CtrlVar.MeshSizeMin);
        end
        
        if CtrlVar.GLmeshing==1
            [coordinates,connectivity]=genmesh2d(Experiment,MeshBoundaryCooWithGLcoo,CtrlVar,edge,face);
        else
            [coordinates,connectivity]=genmesh2d(Experiment,MeshBoundaryCoordinates,CtrlVar,[],[],GmeshBackgroundScalarField);
        end
        
        [Nele,~]=size(connectivity);
        fprintf(CtrlVar.fidlog,'new Nele after rescaling is %-i and CtrlVar.MaxNumberOfElements is %-i  \n',Nele,CtrlVar.MaxNumberOfElements);
    end
    
end

end

