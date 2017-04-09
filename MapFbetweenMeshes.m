function [UserVar,Fnew,BCsNew,GFnew,lnew]=MapFbetweenMeshes(UserVar,CtrlVar,MUAold,MUAnew,Fold,BCsOld,GFold,lold)

         

narginchk(8,8)
nargoutchk(5,5)

isMeshChanged=HasMeshChanged(MUAold,MUAnew);

Fnew=Fold;

if ~isMeshChanged
    GFnew=GFold;
    BCsNew=BCsOld;
    lnew=lold;
    return
end

MUAnew=UpdateMUA(CtrlVar,MUAnew);
lnew=UaLagrangeVariables;

x=MUAnew.coordinates(:,1); y=MUAnew.coordinates(:,2);


if CtrlVar.TimeDependentRun
    
    % if time dependent then surface (s) and bed (b) are defined by mapping old thickness onto
    [UserVar,~,~,Fnew.S,Fnew.B,Fnew.alpha]=GetGeometry(UserVar,CtrlVar,MUAnew,CtrlVar.time,'SB');
    OutsideValue=0;
    Fnew.h=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,x,y,OutsideValue,Fold.h);
    
else
    
    % if a diagnostic step then surface (s) and bed (b), and hence the thickness (h), are defined by the user
    fprintf('Note that as this is not a time-dependent run the ice upper and lower surfaces (s and b) are defined by the user. \n')
    fprintf('When mapping quantities from an old to a new mesh, all geometrical variables (s, b, S, and B) of the new mesh \n')
    fprintf('are therefore obtained through a call to DefineGeometry.m and not through interpolation from the old mesh.\n')
    
    [UserVar,Fnew.s,Fnew.b,Fnew.S,Fnew.B,Fnew.alpha]=GetGeometry(UserVar,CtrlVar,MUAnew,CtrlVar.time,'sbSB');
    
    Fnew.h=Fnew.s-Fnew.b;
    
end

% I define through user input:
%
% S, B, rho, rhow, as, ab, dasdh, dabdh, C, m, AGlen, n, BCs and all sea ice
% parameters.
%
%

[UserVar,Fnew]=GetDensities(UserVar,CtrlVar,MUAnew,Fnew);
[Fnew.b,Fnew.s,Fnew.h,GFnew]=Calc_bs_From_hBS(CtrlVar,MUAnew,Fnew.h,Fnew.S,Fnew.B,Fnew.rho,Fnew.rhow);


[UserVar,Fnew]=GetSlipperyDistribution(UserVar,CtrlVar,MUAnew,Fnew,GFnew);
[UserVar,Fnew]=GetAGlenDistribution(UserVar,CtrlVar,MUAnew,Fnew,GFnew);
[UserVar,Fnew]=GetMassBalance(UserVar,CtrlVar,MUAnew,Fnew,GFnew);

BCsNew=BoundaryConditions;
[UserVar,BCsNew]=GetBoundaryConditions(UserVar,CtrlVar,MUAnew,BCsNew,Fnew,GFnew);

[UserVar,Fnew]=GetSeaIceParameters(UserVar,CtrlVar,MUAnew,BCsNew,Fnew,GFnew);

OutsideValues=[];

[Fnew.ub,Fnew.vb,Fnew.ud,Fnew.vd,Fnew.dhdt,Fnew.dubdt,Fnew.dvbdt,Fnew.duddt,Fnew.dvddt]=...
    MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,x,y,OutsideValues,...
    Fold.ub,Fold.vb,Fold.ud,Fold.vd,Fold.dhdt,Fold.dubdt,Fold.dvbdt,Fold.duddt,Fold.dvddt);


end




