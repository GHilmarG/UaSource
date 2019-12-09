function [UserVar,RunInfo,Fnew,BCsNew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,lold)

         

narginchk(8,8)
nargoutchk(5,5)

RunInfo.MeshAdapt.isChanged=HasMeshChanged(MUAold,MUAnew);

Fnew=Fold;

if ~RunInfo.MeshAdapt.isChanged

    BCsNew=BCsOld;
    lnew=lold;
    return
end

MUAnew=UpdateMUA(CtrlVar,MUAnew);
lnew=UaLagrangeVariables;

Fnew.GF=[] ; % make sure to reset GF if the mesh has changed.  GF can only be calculated once both the new
             % density and the new geometry has been interpolated onto the new mesh. 

x=MUAnew.coordinates(:,1); y=MUAnew.coordinates(:,2);


if CtrlVar.TimeDependentRun
    
    
    if CtrlVar.CurrentRunStepNumber==1  && ~CtrlVar.Restart
        
        fprintf('Note: As this is the first run-step in a time-dependent run: \n')
        fprintf('        When mapping quantities from an old to a new mesh, all geometrical variables (s, b, S, and B) of the new mesh \n')
        fprintf('        are defined through a call to DefineGeometry.m and not through interpolation from the old mesh.\n')
        
        [UserVar,Fnew]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,'sbSB');
        
        %[UserVar,Fnew.s,Fnew.b,Fnew.S,Fnew.B,Fnew.alpha]=GetGeometry(UserVar,CtrlVar,MUAnew,CtrlVar.time,'sbSB');
        %Fnew.h=Fnew.s-Fnew.b;
        
    else
        % if time dependent then surface (s) and bed (b) are defined by mapping old thickness onto
        OutsideValue=[];
        
        if CtrlVar.MapOldToNew.Surface
            
            CtrlVar.Calculate.Geometry="bh-FROM-sBS" ; %    {"bs-FROM-hBS" ; "hb-FROM-sBS" }
            [Fnew.s,Fnew.b]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,x,y,OutsideValue,Fold.s,Fold.b);
            [UserVar,Fnew]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,'SB');
            
            
            CtrlVar.Calculate.Geometry="bs-FROM-hBS" ;
            FTest=Fold; FTest.GF=[]; 
            FTest.h=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,x,y,OutsideValue,FTest.h);
            [UserVar,FTest]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,FTest,'SB');
            
            
            
            FindOrCreateFigure("TestingMapping",[100 100 1000 1000]);
            subplot(2,2,1)
            PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.s) ; title('s : bh-FROM-sBS')
            subplot(2,2,2)
            PlotMeshScalarVariable(CtrlVar,MUAnew,FTest.s); title('s : bs-FROM-hBS')
            subplot(2,2,3)
            PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.h) ; title('h : bh-FROM-sBS')
            subplot(2,2,4)
            PlotMeshScalarVariable(CtrlVar,MUAnew,FTest.h); title('h : bs-FROM-hBS')
            
            
        else
            
            CtrlVar.Calculate.Geometry="bs-FROM-hBS" ;
            Fnew.h=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,x,y,OutsideValue,Fold.h);
            [UserVar,Fnew]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,'SB');
            
        end
        
        
        
        
        % I calculate s and b from h
        % [UserVar,~,~,Fnew.S,Fnew.B,Fnew.alpha]=GetGeometry(UserVar,CtrlVar,MUAnew,CtrlVar.time,'SB');
        % [UserVar,Fnew]=GetDensities(UserVar,CtrlVar,MUAnew,Fnew);
        % [Fnew.b,Fnew.s,Fnew.h,GFnew]=Calc_bs_From_hBS(CtrlVar,MUAnew,Fnew.h,Fnew.S,Fnew.B,Fnew.rho,Fnew.rhow);
    end
    
else
    
    % if a diagnostic step then surface (s) and bed (b), and hence the thickness (h), are defined by the user
    fprintf('Note: As this is not a time-dependent run: \n')
    fprintf('        When mapping quantities from an old to a new mesh, all geometrical variables (s, b, S, and B) of the new mesh \n')
    fprintf('        are defined through a call to DefineGeometry.m and not through interpolation from the old mesh.\n')
    
    [UserVar,Fnew]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,'sbSB');
    
    %[UserVar,Fnew.s,Fnew.b,Fnew.S,Fnew.B,Fnew.alpha]=GetGeometry(UserVar,CtrlVar,MUAnew,CtrlVar.time,'sbSB');
    %Fnew.h=Fnew.s-Fnew.b;
    
end

% I define through user input:
%
% S, B, rho, rhow, as, ab, dasdh, dabdh, C, m, AGlen, n, BCs and all sea ice
% parameters.
%
%

%[UserVar,Fnew]=GetDensities(UserVar,CtrlVar,MUAnew,Fnew);
%[Fnew.b,Fnew.s,Fnew.h,GFnew]=Calc_bs_From_hBS(CtrlVar,MUAnew,Fnew.h,Fnew.S,Fnew.B,Fnew.rho,Fnew.rhow);


[UserVar,Fnew]=GetSlipperyDistribution(UserVar,CtrlVar,MUAnew,Fnew);
[UserVar,Fnew]=GetAGlenDistribution(UserVar,CtrlVar,MUAnew,Fnew);
[UserVar,Fnew]=GetMassBalance(UserVar,CtrlVar,MUAnew,Fnew);

BCsNew=BoundaryConditions;
[UserVar,BCsNew]=GetBoundaryConditions(UserVar,CtrlVar,MUAnew,BCsNew,Fnew);

[UserVar,Fnew]=GetSeaIceParameters(UserVar,CtrlVar,MUAnew,BCsNew,Fnew);

OutsideValues=[];

%%


[Fnew.ub,Fnew.vb,Fnew.ud,Fnew.vd,Fnew.dhdt,Fnew.dubdt,Fnew.dvbdt,Fnew.duddt,Fnew.dvddt]=...
    MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,x,y,OutsideValues,...
    Fold.ub,Fold.vb,Fold.ud,Fold.vd,Fold.dhdt,Fold.dubdt,Fold.dvbdt,Fold.duddt,Fold.dvddt);



%%
end




