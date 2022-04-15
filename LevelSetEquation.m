function [UserVar,RunInfo,LSF,LSFMask,LSFnodes,l,LSFqx,LSFqy]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)
%%
%
%
%    df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f0)
%
%    df/dt + (u-cx) df/dx + (v-cy) df/dy - div (kappa grad f) = 0
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

if any(isnan(F0.c))
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


    % DistEle=pdist2([xc(:) yc(:)],[MUA.xEle MUA.yEle],'euclidean','Smallest',1) ;
    % DistEle=DistEle(:) ;  % note, this is a element-valued distance function
    F0.x=MUA.coordinates(:,1); F0.y=MUA.coordinates(:,2);
    DistNod=pdist2([xc(:) yc(:)],[F0.x F0.y],'euclidean','Smallest',1) ;
    DistNod=DistNod(:) ;  % note, this is a element-valued distance function
    DistEle=Nodes2EleMean(MUA.connectivity,DistNod) ;

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

    F1.LSF=F1.LSF(kk) ;
    F1.ub=F1.ub(kk);
    F1.vb=F1.vb(kk);



    % additonal variables for sliding law evaluation at int point
    F1.h=F1.h(kk) ;      F0.h=F0.h(kk) ;
    F1.s=F1.s(kk) ;      F0.s=F0.s(kk) ;
    F1.b=F1.b(kk) ;      F0.b=F0.b(kk) ;
    F1.S=F1.S(kk) ;      F0.S=F0.S(kk) ;
    F1.rho=F1.rho(kk) ;  F0.rho=F0.rho(kk) ;

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


[UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationSolver(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);


%
% FBCs=FindOrCreateFigure("BCs LSF") ; clf(FBCs) ;
% PlotBoundaryConditions(CtrlVar,MUA,BCs) ;

if CtrlVar.LevelSetMethodSolveOnAStrip


    FindOrCreateFigure("LSF solve mesh ") ; PlotMuaMesh(CtrlVar,MUA) ;

    Flsf=FindOrCreateFigure("LSF on strip ") ; clf(Flsf) ;
    PlotMeshScalarVariable(CtrlVar,MUA,LSF/1000) ;
    hold on ; PlotCalvingFronts(CtrlVar,MUA,LSF,"r",LineWidth=2 ) ;
    title("$\varphi$"+sprintf(" at t=%2.2f",F1.time),Interpreter="latex")


    LSFcopy(~LSFnodes)=DistNod(~LSFnodes).*sign(LSFcopy(~LSFnodes)) ;  % over old nodes at save distance from zero-level, fill in using already calculated signed distance
    % as these values will never impact the calculation of the zero level

    %LSFcopy(LSFnodes)=LSF;

    LSFcopy(LSFnodes)=LSF(ll(LSFnodes));

   % LSFcopy(~isnan(ll))=LSF(ll(~isnan(ll)));



    Flsf=FindOrCreateFigure("LSF on original mesh ") ; clf(Flsf) ;
    PlotMeshScalarVariable(CtrlVar,MUAonInput,LSFcopy/1000) ;
    hold on ; PlotCalvingFronts(CtrlVar,MUAonInput,LSFcopy,"r",LineWidth=2 ) ;
    title("$\varphi$"+sprintf(" at t=%2.2f",F1.time),Interpreter="latex")

    LSF=LSFcopy;

end

LSFMask=CalcMeshMask(CtrlVar,MUAonInput,LSF,0);  % If I solved the LSF on a strip, this will not be the correct mask over the full MUA
% unless I use the original MUA




% F1.LSF=LSF;
% F1.LSFMask=Mask;  % If I solved the LSF on a strip, this will not be the correct mask over the full MUA
% F1.LSFqx=LSFqx;
% F1.LSFqy=LSFqy;
% F1.LSFnodes=LSFnodes;



end



