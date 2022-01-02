function [UserVar,RunInfo,LSF,Mask,l,LSFqx,LSFqy]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)
%%
%
%
%    df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f0)
%
%    df/dt + (u-cx) df/dx + (v-cy) df/dy - div (kappa grad f) = 0
%
%

narginchk(7,8)
nargoutchk(7,7)



if ~CtrlVar.LevelSetMethod

    LSF=[];
    l=[];
    LSFqx=[];
    LSFqy=[];
    return
end

if any(isnan(F0.c))
    fprintf("Level set is not evolved because calving rate (c) contains nan. \n")
    LSF=F1.LSF;
    Mask=[];
    l=[];
    LSFqx=[];
    LSFqy=[];
    return
end


if nargin<8
    l=[];
end




if CtrlVar.LevelSetMethodSolveOnAStrip

 

    CtrlVar.LineUpGLs=false ; Threshold=0 ;

    [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F0.LSF,Threshold);


    % DistEle=pdist2([xc(:) yc(:)],[MUA.xEle MUA.yEle],'euclidean','Smallest',1) ;
    % DistEle=DistEle(:) ;  % note, this is a element-valued distance function

    DistNod=pdist2([xc(:) yc(:)],[F0.x F0.y],'euclidean','Smallest',1) ;
    DistNod=DistNod(:) ;  % note, this is a element-valued distance function
    DistEle=Nodes2EleMean(MUA.connectivity,DistNod) ;


  
    ElementsToBeDeactivated=DistEle>CtrlVar.LevelSetMethodStripWidth;

    [MUA,K]=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated)  ;

    LSFcopy=F0.LSF ; % make a copy of LSF 
    F0.LSF=F0.LSF(K) ;
    F0.ub=F0.ub(K);
    F0.vb=F0.vb(K);

    F1.LSF=F1.LSF(K) ;
    F1.ub=F1.ub(K);
    F1.vb=F1.vb(K);

    if ~isempty(F0.c)
        F0.c=F0.c(K);
    end

    if  ~isempty(F1.c)
        F1.c=F1.c(K);
    end



end


[UserVar,RunInfo,LSF,Mask,l,LSFqx,LSFqy]=LevelSetEquationSolver(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);

FindOrCreateFigure("LSF solve mesh ") ; PlotMuaMesh(CtrlVar,MUA) ;
Flsf=FindOrCreateFigure("LSF solve LSF ") ; clf(Flsf) ;
PlotMeshScalarVariable(CtrlVar,MUA,LSF/1000) ;
hold on ; PlotCalvingFronts(CtrlVar,MUA,LSF,"r",LineWidth=2 ) ; 
title("$\varphi$"+sprintf(" at t=%2.1f",F1.time),Interpreter="latex")

if CtrlVar.LevelSetMethodSolveOnAStrip
    
    % 
    io=LSFcopy < 0 ;
    LSFcopy(io)=-DistNod(io)  ;   % were not evolved, just use the signed distance that has already been calculated
    LSFcopy(K)=LSF;               % and then update the LSF field where it was evolved
    LSF=LSFcopy;   


end





end



