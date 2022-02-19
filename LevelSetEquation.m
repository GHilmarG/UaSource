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
    LSFnodes=~isnan(ll) ; % a logical list with the size MUA.Nnodes, true for nodes in MUA over which the level-set is evolved



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
    F1.h=F1.h(kk) ;
    F0.h=F0.h(kk) ;

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

% FindOrCreateFigure("LSF solve mesh ") ; PlotMuaMesh(CtrlVar,MUA) ;
%
% Flsf=FindOrCreateFigure("LSF solve LSF ") ; clf(Flsf) ;
% PlotMeshScalarVariable(CtrlVar,MUA,LSF/1000) ;
% hold on ; PlotCalvingFronts(CtrlVar,MUA,LSF,"r",LineWidth=2 ) ;
% title("$\varphi$"+sprintf(" at t=%2.2f",F1.time),Interpreter="latex")
%
% FBCs=FindOrCreateFigure("BCs LSF") ; clf(FBCs) ;
% PlotBoundaryConditions(CtrlVar,MUA,BCs) ;

if CtrlVar.LevelSetMethodSolveOnAStrip
    %     %
    %     io=LSFcopy < 0 ;
    %     LSFcopy(io)=-DistNod(io)  ;   % were not evolved, just use the signed distance that has already been calculated
    %     LSFcopy(kk)=LSF;               % and then update the LSF field where it was evolved

    % new
    io=isnan(ll);                                 % where I have only the old (before deactivation) nodes,
    LSFcopy(io)=DistNod(io).*sign(LSFcopy(io)) ;  % over old nodes at save distance from zero-level, fill in using already calculated signed distance
    % as these values will never impact the calculation of the zero level
    LSFcopy(~io)=LSF;                             % over the strip around zero level, use the calculated values based on solving the LSF equation

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



