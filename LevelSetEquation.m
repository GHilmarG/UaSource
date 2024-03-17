function [UserVar,RunInfo,LSF,LSFMask,LSFnodes,l,LSFqx,LSFqy]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)
%%
%
% 
% 
% $$\partial_t f +  \mathbf{v} \cdot \nabla f  - \nabla \cdot (\kappa \nabla f) = c \, \|(\nabla f)\|$$
% 
%
%  df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f)
%
%    df/dt + (u-cx) df/dx + (v-cy) df/dy - div (kappa grad f) = 0
%
%%

%
% Note: Here, F0 and F1 have both been calculated in a uvh solve. 
%       The F1 contains the solve at time F1.time, which is now the `current time', ie the time at which uvh has been solved
%       
%       

narginchk(7,8)
nargoutchk(8,8)

LSF=[];
l=[];
LSFqx=[];
LSFqy=[];
LSFMask=[];
LSFnodes=[];

if ~CtrlVar.LevelSetMethod
    return
end

if any(isnan(F0.c)) || isempty(F0.c)
    fprintf("Level set is not evolved because calving rate (c) contains nan. \n")
    LSF=F1.LSF;
    return
end


if nargin<8
    l=[];
end

MUAonInput=MUA;  % I need this in case I do the solution on a subset


if CtrlVar.LevelSetMethodSolveOnAStrip



    CtrlVar.LineUpGLs=false ; Threshold=0 ;

    [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F0.LSF,Threshold);
    F0.x=MUA.coordinates(:,1); F0.y=MUA.coordinates(:,2);
    DistNod=pdist2([xc(:) yc(:)],[F0.x F0.y],'euclidean','Smallest',1) ;
    DistNod=DistNod(:) ;  
    DistEle=Nodes2EleMean(MUA.connectivity,DistNod) ; % note, this is now an element-valued distance function

    if isnan(CtrlVar.LevelSetMethodStripWidth)

        fprintf("The variable CtrlVar.LevelSetMethodStripWidth needs to be defined.\n")
        error("LevelSetEquation:ParameterNotDefined","The variable CtrlVar.LevelSetMethodStripWidth needs to be defined.")

    end


    ElementsToBeDeactivated=DistEle>CtrlVar.LevelSetMethodStripWidth;

    [MUA,kk,ll]=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated)  ;
    LSFnodes=~isnan(ll) ; % a logical list with the size MUA.Nnodes, true for nodes in MUA (i.e. the MUA here given as an input)
    % over which the level-set is evolved

    % Thist is a bit of a lazy approach because I know which nodes were deleted and the new nodal numbers
    % so it would be possibly to figure out how to map the BCs from the old to new.
    % Also, some of the issues are just related to the uvh boundary conditions that are not relevant when solving for \varphi

    % [UserVar,BCs]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,F0) ;
    BCs=BoundaryConditions;  % resetting

    LSFcopy=F0.LSF ; % make a copy of LSF
    F0.LSF=F0.LSF(kk) ;
    F0.ub=F0.ub(kk);
    F0.vb=F0.vb(kk);

    F0.LSFMask.NodesIn=F0.LSFMask.NodesIn(kk);
    F0.LSFMask.NodesOn=F0.LSFMask.NodesOn(kk);
    F0.LSFMask.NodesOut=F0.LSFMask.NodesOut(kk);

    F1.LSF=F1.LSF(kk) ;
    F1.ub=F1.ub(kk);
    F1.vb=F1.vb(kk);

    % additonal variables for sliding law evaluation at int point
    F1.x=F1.x(kk) ;      F0.x=F0.x(kk) ;
    F1.y=F1.y(kk) ;      F0.y=F0.y(kk) ;
    F1.h=F1.h(kk) ;      F0.h=F0.h(kk) ;
    F1.s=F1.s(kk) ;      F0.s=F0.s(kk) ;
    F1.b=F1.b(kk) ;      F0.b=F0.b(kk) ;
    F1.S=F1.S(kk) ;      F0.S=F0.S(kk) ;
    F1.B=F1.B(kk) ;      F0.B=F0.B(kk) ;

    F1.C=F1.C(kk) ;      F0.C=F0.C(kk) ;
    F1.AGlen=F1.AGlen(kk) ;      F0.AGlen=F0.AGlen(kk) ;
    F1.GF.node=F1.GF.node(kk) ;      F0.GF.node=F0.GF.node(kk) ;


    % To do, set all other values to empty to make sure they are not updated

    if ~isempty(F0.c)
        F0.c=F0.c(kk);
    end

    if  ~isempty(F1.c)
        F1.c=F1.c(kk);
    end



    % How to update BCs? I need to have a mapping from old-to-new nodal numbers

else
    LSFnodes=[] ;
end


%% TestIng: Calculating  various potential calving-law related quantities ahead of a call to the level-set equation solver
if CtrlVar.LevelSetMethodTest 

    [F0.exx,F0.eyy,F0.exy]=CalcNodalStrainRates(MUA,F0.ub,F0.vb);
    [F1.exx,F1.eyy,F1.exy]=CalcNodalStrainRates(MUA,F1.ub,F1.vb);

    PSR=CalcPrincipalValuesOfSymmetricalTwoByTwoMatrices(F0.exx,F0.exy,F0.eyy); % Principal Strain Rates
    I1=PSR(:,1)<0 ;  PSR(I1,1)=0;
    I2=PSR(:,2)<0 ;  PSR(I2,2)=0;

    K=1;
    cEigenCalving=K*PSR(:,1).*PSR(:,2);


    PSR(F0.LSFMask.NodesOut,:)=NaN;
    FindOrCreateFigure("P1") ; PlotMeshScalarVariable(CtrlVar,MUA,PSR(:,1)) ;
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
    hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w') ;
    hold on ; PlotCalvingFronts(CtrlVar,MUA,F0,'r');

    FindOrCreateFigure("P2") ; PlotMeshScalarVariable(CtrlVar,MUA,PSR(:,2)) ;
    hold on ; PlotCalvingFronts(CtrlVar,MUA,F0,'r');



    FindOrCreateFigure("Eigen Calving") ; PlotMeshScalarVariable(CtrlVar,MUA,cEigenCalving) ;

    scale=1 ; FindOrCreateFigure("strain rates F0"); PlotTensor(F0.x/1000,F0.y/1000,F0.exx,F0.exy,F0.eyy,scale) ;  axis equal
    hold on ; PlotCalvingFronts(CtrlVar,MUA,F0,'r');
    scale=1 ; FindOrCreateFigure("strain rates F1"); PlotTensor(F1.x/1000,F1.y/1000,F1.exx,F1.exy,F1.eyy,scale) ;  axis equal
    hold on ; PlotCalvingFronts(CtrlVar,MUA,F0,'r');

end

%%

[UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationSolver(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);


if CtrlVar.LevelSetMethodSolveOnAStrip


    if CtrlVar.LevelSetInfoLevel>=10 && CtrlVar.doplots


        Flsf=FindOrCreateFigure("LSF on strip ") ; clf(Flsf) ;
        PlotMeshScalarVariable(CtrlVar,MUA,LSF/1000) ;
        hold on ; PlotCalvingFronts(CtrlVar,MUA,LSF,"r",LineWidth=2 ) ;
        title("$\varphi$"+sprintf(" at t=%2.2f",F1.time),Interpreter="latex")

    end


    LSFcopy(~LSFnodes)=DistNod(~LSFnodes).*sign(LSFcopy(~LSFnodes)) ;  
    % over old nodes at save distance from zero-level, fill in using already calculated signed distance
    % as these values will never impact the calculation of the zero level
    LSFcopy(LSFnodes)=LSF(ll(LSFnodes));

    % LSFcopy(~isnan(ll))=LSF(ll(~isnan(ll)));


    if CtrlVar.LevelSetInfoLevel>=10 && CtrlVar.doplots
        Flsf=FindOrCreateFigure("LSF on original mesh ") ; clf(Flsf) ;
        PlotMeshScalarVariable(CtrlVar,MUAonInput,LSFcopy/1000) ;
        hold on ; PlotCalvingFronts(CtrlVar,MUAonInput,LSFcopy,"r",LineWidth=2 ) ;
        title("$\varphi$"+sprintf(" at t=%2.2f",F1.time),Interpreter="latex")
    end

    LSF=LSFcopy;

end

LSFMask=CalcMeshMask(CtrlVar,MUAonInput,LSF,0);  % If I solved the LSF on a strip, this will not be the correct mask over the full MUA
                                                 % unless I use the original MUA

end



