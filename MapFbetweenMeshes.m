function [UserVar,RunInfo,Fnew,BCsNew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,lold,OutsideValue)

% �a


narginchk(8,9)
nargoutchk(3,5)

nouts=nargout;

if nargin<9 || isempty(OutsideValue)
    OutsideValue.h=CtrlVar.ThickMin;
    OutsideValue.b=0;
    OutsideValue.s=CtrlVar.ThickMin;
end

RunInfo.MeshAdapt.isChanged=HasMeshChanged(MUAold,MUAnew);

% Fnew=Fold;



if ~RunInfo.MeshAdapt.isChanged
    Fnew=Fold;
    BCsNew=BCsOld;
    lnew=lold;

    RunInfo.Mapping.nNewNodes=MUAnew.Nnodes;
    RunInfo.Mapping.nOldNodes=MUAold.Nnodes;
    RunInfo.Mapping.nIdenticalNodes=MUAnew.Nnodes;
    RunInfo.Mapping.nNotIdenticalNodes=0;
    RunInfo.Mapping.nNotIdenticalNodesOutside=0;
    RunInfo.Mapping.nNotIdenticalNodesInside=0 ;
    RunInfo.Mapping.nNotIdenticalInside=0;

    return
end

MUAnew=UpdateMUA(CtrlVar,MUAnew);
lnew=UaLagrangeVariables;

Fnew=UaFields;
Fnew.time=Fold.time;
Fnew.GF=[] ; % make sure to reset GF if the mesh has changed.  GF can only be calculated once both the new
% density and the new geometry has been interpolated onto the new mesh.

Fnew.x=MUAnew.coordinates(:,1); Fnew.y=MUAnew.coordinates(:,2);
Fold.x=MUAold.coordinates(:,1); Fold.y=MUAold.coordinates(:,2);


if CtrlVar.TimeDependentRun


    if CtrlVar.CurrentRunStepNumber==1  && ~CtrlVar.Restart

        fprintf('Note: As this is the first run-step in a time-dependent run: \n')
        fprintf('        When mapping quantities from an old to a new mesh, all geometrical variables (s, b, S, and B) of the new mesh \n')
        fprintf('        are defined through a call to DefineGeometry.m and not through interpolation from the old mesh.\n')

        [UserVar,Fnew]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,'-s-b-S-B-rho-rhow-g-');

        if  CtrlVar.LevelSetMethod

            Fnew.LSF=[] ; % force a re-calculation of the level-set
            fprintf('        The level set is also now reset to empty, and then defined by a call to DefineCalving.\n')
            [UserVar,Fnew]=GetCalving(UserVar,CtrlVar,MUAnew,Fnew,[]);

            fprintf("Setting ice thicknesses downstream of calving fronts to the minimum prescribed value of %f .\n",CtrlVar.LevelSetMinIceThickness)
            Fnew.h(Fnew.LSFMask.NodesOut)=CtrlVar.LevelSetMinIceThickness;
            [Fnew.b,Fnew.s,Fnew.h,Fnew.GF]=Calc_bs_From_hBS(CtrlVar,MUAnew,Fnew.h,Fnew.S,Fnew.B,Fnew.rho,Fnew.rhow);

        end

    else  % this is time-dependent and not start of a new run
          % here s and b must be interpolated and LSF if using the level-set method


        switch CtrlVar.MapOldToNew.Transient.Geometry
            
            % This is the default approach in a transient simulation.  Here only s is interpolated onto the new mesh. B and S are then
            % later obtained through a call to GetGeometryAndDensities and in that call, b and h are calculated from s, S and B.
            %
            % This approach does not conserve h!
            %
            % This is done quite deliberately because the surface s is genearly very smooth, whereas b and h are not. If, for example,
            % the new mesh has grounded nodes located over a local topographic feature in the bed, and h was interpolated from old to new
            % mesh, the new bedrock topographic feature would impact the surface s=B+h and produce unrealistic surface topography.

            case "bh-FROM-sBS"


                if ~isfield(OutsideValue,'s')
                    OutsideValue.s=NaN;
                end
                if ~isfield(OutsideValue,'b')
                    OutsideValue.b=NaN;
                end

                [RunInfo,Fnew.s]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValue.s,Fold.s);


            case "bs-FROM-hBS"


                if ~isfield(OutsideValue,'h')
                    OutsideValue.h=NaN;
                end

                [RunInfo,Fnew.h]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValue.h,Fold.h);

        end


        CtrlVar.Calculate.Geometry=CtrlVar.MapOldToNew.Transient.Geometry ;
        % Here within GetGeometryAndDensities, here b and h are calculated from s, S and B
        [UserVar,Fnew]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,"-S-B-rho-");
        % Important to do this ahead of call the GetCalving in case the user wants to define LSF in terms of GF (as points out by
        % Sainan Sun on 2nd August, 2023)


        if  CtrlVar.LevelSetMethod

            if CtrlVar.LevelSetEvolution=="-Prescribed-"

                fprintf("MapFbetweenMeshes: LevelSetEvolution is prescribed, so when mapping onto a new mesh, the levelset is defined through a call to DefineCalving.m \n")
                BCsNew=[] ; % BCs have yet to be defined
                Fnew.LSF=[] ;
                [UserVar,Fnew]=GetCalving(UserVar,CtrlVar,MUAnew,Fnew,BCsNew) ;

            else


                if ~isfield(OutsideValue,'LSF')
                    OutsideValue.LSF=NaN;
                end

                [RunInfo,Fnew.LSF]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValue.LSF,Fold.LSF);

            end

        end




        if CtrlVar.MapOldToNew.Test

            % calculate geometry using both options

            OutsideValue=[];
            [RunInfo,Fnew.s,Fnew.b]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValue,Fold.s,Fold.b);
            CtrlVar.Calculate.Geometry="bh-FROM-sBS" ;

            [UserVar,Fnew]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,"-S-B-rho-");

            FTest=Fold; FTest.GF=[];
            [RunInfo,FTest.h]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValue,Fold.h);
            CtrlVar.Calculate.Geometry="bs-FROM-hBS" ;
            [UserVar,FTest]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,FTest,"-S-B-rho-");

            % now plot

            FindOrCreateFigure("TestingMapping");
            hold off

            subplot(3,4,1) ; PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.s) ; title('s : bh-FROM-sBS')
            hold on ; PlotGroundingLines(CtrlVar,MUAnew,Fnew.GF);
            subplot(3,4,2) ; PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.b) ; title('b : bh-FROM-sBS')
            subplot(3,4,3) ; PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.B) ; title('B : bh-FROM-sBS')
            subplot(3,4,4) ; PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.h) ; title('h : bh-FROM-sBS')


            subplot(3,4,5) ; PlotMeshScalarVariable(CtrlVar,MUAnew,FTest.s); title('s : bs-FROM-hBS')
            hold on ;        PlotGroundingLines(CtrlVar,MUAnew,FTest.GF);
            subplot(3,4,6) ; PlotMeshScalarVariable(CtrlVar,MUAnew,FTest.b); title('b : bs-FROM-hBS')
            subplot(3,4,7) ; PlotMeshScalarVariable(CtrlVar,MUAnew,FTest.B); title('B : bs-FROM-hBS')
            subplot(3,4,8) ; PlotMeshScalarVariable(CtrlVar,MUAnew,FTest.h); title('h : bs-FROM-hBS')

            subplot(3,4,9) ; PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.s-FTest.s); title('ds')
            subplot(3,4,10) ; PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.b-FTest.b); title('db')
            subplot(3,4,11) ; PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.B-FTest.B); title('dB')
            subplot(3,4,12) ; PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.h-FTest.h); title('dh')

        end

    end

else

    % if a diagnostic step then surface (s) and bed (b), and hence the thickness (h), are defined by the user
    fprintf('Note: As this is not a time-dependent run: \n')
    fprintf('        When mapping quantities from an old to a new mesh, all geometrical variables (s, b, S, and B) of the new mesh \n')
    fprintf('        are defined through a call to DefineGeometry and not through interpolation from the old mesh.\n')


    [UserVar,Fnew]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,'-s-b-S-B-rho-rhow-g-');

    if  CtrlVar.LevelSetMethod
        Fnew.LSF=[] ; % force a re-calculation of the level-set
        fprintf('        The level set is also now reset to empty,and then defined by a call to DefineCalving.\n')
        [UserVar,Fnew]=GetCalving(UserVar,CtrlVar,MUAnew,Fnew,[]);

        fprintf("Setting ice thicknesses downstream of calving fronts to the minimum prescribed value of %f .\n",CtrlVar.LevelSetMinIceThickness)
        Fnew.h(Fnew.LSFMask.NodesOut)=CtrlVar.LevelSetMinIceThickness;
        [Fnew.b,Fnew.s,Fnew.h,Fnew.GF]=Calc_bs_From_hBS(CtrlVar,MUAnew,Fnew.h,Fnew.S,Fnew.B,Fnew.rho,Fnew.rhow);

    end


end

Fnew.x=MUAnew.coordinates(:,1) ;  Fnew.y=MUAnew.coordinates(:,2) ;

[UserVar,Fnew]=GetSlipperyDistribution(UserVar,CtrlVar,MUAnew,Fnew);


[UserVar,Fnew]=GetAGlenDistribution(UserVar,CtrlVar,MUAnew,Fnew);
[UserVar,Fnew]=GetMassBalance(UserVar,CtrlVar,MUAnew,Fnew);


if nouts>=4
    BCsNew=BoundaryConditions;
    [UserVar,BCsNew]=GetBoundaryConditions(UserVar,CtrlVar,MUAnew,BCsNew,Fnew);
else
    BCsNew=[];
end

% [UserVar,Fnew]=GetCalving(UserVar,CtrlVar,MUAnew,Fnew,BCsNew) ;


[UserVar,Fnew]=GetSeaIceParameters(UserVar,CtrlVar,MUAnew,BCsNew,Fnew);

%%


if ~isfield(OutsideValue,'ub')
    OutsideValue.ub=NaN;
end
if ~isfield(OutsideValue,'vb')
    OutsideValue.vb=NaN;
end

if ~isfield(OutsideValue,'ud')
    OutsideValue.ud=NaN;
end

if ~isfield(OutsideValue,'vd')
    OutsideValue.vd=NaN;
end

if ~isfield(OutsideValue,'dhdt')
    OutsideValue.dhdt=NaN;
end


if ~isfield(OutsideValue,'dubdt')
    OutsideValue.dubdt=NaN;
end

if ~isfield(OutsideValue,'dvbdt')
    OutsideValue.dvbdt=NaN;
end

if ~isfield(OutsideValue,'duddt')
    OutsideValue.duddt=NaN;
end

if ~isfield(OutsideValue,'dvddt')
    OutsideValue.dvddt=NaN;
end

if ~isfield(OutsideValue,'LSF')
    OutsideValue.LSF=NaN;
end


switch lower(CtrlVar.FlowApproximation)

    case "sstream"

        % if ub vb has been calculated as a part of the remeshing, I'm loosing that information
        % here. The problem is that ub vb will not, in general, have been calculated on either
        % MUAnew or MUAold unless only one adapt iteration was performed.
        %
        [RunInfo,Fnew.ub,Fnew.vb,Fnew.ud,Fnew.vd,Fnew.dhdt,Fnew.dubdt,Fnew.dvbdt]=...
            MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,...
            [OutsideValue.ub,OutsideValue.vb,OutsideValue.ud,OutsideValue.ud,OutsideValue.dhdt,OutsideValue.dubdt,OutsideValue.dvbdt],...
            Fold.ub,Fold.vb,Fold.ud,Fold.vd,Fold.dhdt,Fold.dubdt,Fold.dvbdt) ;


        Fnew.duddt=zeros(MUAnew.Nnodes,1);
        Fnew.dvddt=zeros(MUAnew.Nnodes,1);

    case "ssheet"

        [RunInfo,Fnew.ub,Fnew.vb,Fnew.ud,Fnew.vd,Fnew.dhdt,Fnew.duddt,Fnew.dvddt,Fnew.LSF]=...
            MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,...
            [OutsideValue.ub,OutsideValue.vb,OutsideValue.ud,OutsideValue.ud,OutsideValue.dhdt,OutsideValue.duddt,OutsideValue.dvddt,OutsideValue.LSF],...
            Fold.ub,Fold.vb,Fold.ud,Fold.vd,Fold.dhdt,Fold.duddt,Fold.dvddt,Fold.LSF) ;

        Fnew.dubdt=zeros(MUAnew.Nnodes,1);
        Fnew.dvbdt=zeros(MUAnew.Nnodes,1);

    case "uvhprescribed"

        % This will update velocities and thickness
        [UserVar,RunInfo,Fnew]=uvhPrescibed(UserVar,RunInfo,CtrlVar,MUAnew,Fnew,Fnew,lnew,BCsNew);

        Fnew.dubdt=zeros(MUAnew.Nnodes,1); Fnew.dvbdt=zeros(MUAnew.Nnodes,1);
        Fnew.duddt=zeros(MUAnew.Nnodes,1); Fnew.dvddt=zeros(MUAnew.Nnodes,1);



    otherwise

        error("MapFbetweenMeshes:CaseNotFound","case not found")
end






%%
end




