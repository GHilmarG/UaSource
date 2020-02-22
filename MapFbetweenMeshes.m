function [UserVar,RunInfo,Fnew,BCsNew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,lold,OutsideValue)

         

narginchk(8,9)
nargoutchk(5,5)

if nargin<9 || isempty(OutsideValue)
    OutsideValue.h=CtrlVar.ThickMin;
    OutsideValue.b=0;
    OutsideValue.s=CtrlVar.ThickMin;
end

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
        
        
        switch CtrlVar.MapOldToNew.Transient.Geometry
            
            case "bh-FROM-sBS"
                
                
                if ~isfield(OutsideValue,'s')
                    OutsideValue.s=NaN;
                end
                if ~isfield(OutsideValue,'b')
                    OutsideValue.b=NaN;
                end
                
                [RunInfo,Fnew.s,Fnew.b]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,[OutsideValue.s OutsideValue.b],Fold.s,Fold.b);
                % I don't really need to map b from old to new, but I use this later as an initial
                % guess.
                
                
            case "bs-FROM-hBS"
                
                       
                if ~isfield(OutsideValue,'h')
                    OutsideValue.h=NaN;
                end
                
                [RunInfo,Fnew.h]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValue.h,Fold.h);
                
        end
        
        CtrlVar.Calculate.Geometry=CtrlVar.MapOldToNew.Transient.Geometry ;
        [UserVar,Fnew]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,'SB');
        
        
        if CtrlVar.MapOldToNew.Test
            
            % calculate geometry using both options
            
            OutsideValue=[];
            [RunInfo,Fnew.s,Fnew.b]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValue,Fold.s,Fold.b);
            CtrlVar.Calculate.Geometry="bh-FROM-sBS" ;
            [UserVar,Fnew]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,'SB');
            
            FTest=Fold; FTest.GF=[];
            [RunInfo,FTest.h]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValue,Fold.h);
            CtrlVar.Calculate.Geometry="bs-FROM-hBS" ;
            [UserVar,FTest]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,FTest,'SB');
            
            % now plot
            
            FindOrCreateFigure("TestingMapping");
            hold off
            
            subplot(3,4,1) ; PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.s) ; title('s : bh-FROM-sBS')
            hold on ; [xGL,yGL]=PlotGroundingLines(CtrlVar,MUAnew,Fnew.GF);
            subplot(3,4,2) ; PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.b) ; title('b : bh-FROM-sBS')
            subplot(3,4,3) ; PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.B) ; title('B : bh-FROM-sBS')
            subplot(3,4,4) ; PlotMeshScalarVariable(CtrlVar,MUAnew,Fnew.h) ; title('h : bh-FROM-sBS')
            
            
            subplot(3,4,5) ; PlotMeshScalarVariable(CtrlVar,MUAnew,FTest.s); title('s : bs-FROM-hBS')
            hold on ; [xGL,yGL]=PlotGroundingLines(CtrlVar,MUAnew,FTest.GF);
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


switch CtrlVar.FlowApproximation
    
    case "SSTREAM"
        
        
        [RunInfo,Fnew.ub,Fnew.vb,Fnew.ud,Fnew.vd,Fnew.dhdt,Fnew.dubdt,Fnew.dvbdt]=...
            MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,...
            [OutsideValue.ub,OutsideValue.vb,OutsideValue.ud,OutsideValue.ud,OutsideValue.dhdt,OutsideValue.dubdt,OutsideValue.dvbdt],...
            Fold.ub,Fold.vb,Fold.ud,Fold.vd,Fold.dhdt,Fold.dubdt,Fold.dvbdt) ;
        
        Fnew.duddt=zeros(MUAnew.Nnodes,1);
        Fnew.dvddt=zeros(MUAnew.Nnodes,1);
        
    case "SSHEET"
        
        [RunInfo,Fnew.ub,Fnew.vb,Fnew.ud,Fnew.vd,Fnew.dhdt,Fnew.duddt,Fnew.dvddt]=...
            MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,...
            [OutsideValue.ub,OutsideValue.vb,OutsideValue.ud,OutsideValue.ud,OutsideValue.dhdt,OutsideValue.duddt,OutsideValue.dvddt],...
            Fold.ub,Fold.vb,Fold.ud,Fold.vd,Fold.dhdt,Fold.duddt,Fold,dvddt) ;
        
        Fnew.dubdt=zeros(MUAnew.Nnodes,1);
        Fnew.dvbdt=zeros(MUAnew.Nnodes,1);
        
        
end

%%
end




