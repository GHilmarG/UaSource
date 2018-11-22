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
        OutsideValue=0;
        Fnew.h=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,x,y,OutsideValue,Fold.h);
        [UserVar,Fnew]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,'SB'); 
        
        
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

[Fnew.ub,Fnew.vb,Fnew.ud,Fnew.vd,Fnew.dhdt,Fnew.dubdt,Fnew.dvbdt,Fnew.duddt,Fnew.dvddt]=...
    MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,x,y,OutsideValues,...
    Fold.ub,Fold.vb,Fold.ud,Fold.vd,Fold.dhdt,Fold.dubdt,Fold.dvbdt,Fold.duddt,Fold.dvddt);


end




