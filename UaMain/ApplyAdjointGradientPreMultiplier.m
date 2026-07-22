





function varargout=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,BCsAdjoint,varargin)

narginchk(4,inf)


if ~isa(BCsAdjoint,"BoundaryConditions")
    error("BCs required")
end

varargout=varargin;


%% Using BCs when calculating/projecting the gradient
%
% This is not finalized, and is here as a test. The uv BCs are already used in the adjoint calculation, the question is if the
% h-BCs should also be used in combination with the use of dh/dt in the data. Something to think more about.
UseBCs=true;
UseBCs=false;

%%



if CtrlVar.Inverse.AdjointGradientPreMultiplier=="M"
    if isa(MUA.dM,"decomposition") && ~UseBCs  % No need for the decomposition object in combination with the BCs, because the KKT system does not use
        % the dM.
        P=MUA.dM/MUA.Area ;
    else
        P=MUA.M/MUA.Area ;
    end



elseif CtrlVar.Inverse.AdjointGradientPreMultiplier=="D"

    l=1e-10;
    P=MUA.Dxx+MUA.Dyy+l*MUA.M;

elseif CtrlVar.Inverse.AdjointGradientPreMultiplier=="H1"

    CtrlVar.Inverse.AdjointGradientPreMultiplier="H1"; 
    ga=CtrlVar.Inverse.PreMultiplier.H1.ga;
    gs=CtrlVar.Inverse.PreMultiplier.H1.gs;

    
    P=gs*(MUA.Dxx+MUA.Dyy)+ga*MUA.M; 

   


elseif CtrlVar.Inverse.AdjointGradientPreMultiplier=="I" ...
        || CtrlVar.Inverse.AdjointGradientPreMultiplier=="Hanalytical"

    varargout=varargin;

    return
end



for k=1:numel(varargin)

    if ~isempty(varargin{k})

        if UseBCs


            [hL,hRhs]=createLh(MUA.Nnodes,BCsAdjoint.hFixedNode,BCsAdjoint.hFixedValue,BCsAdjoint.hTiedNodeA,BCsAdjoint.hTiedNodeB);
            f=varargin{k};
            x0=[]; y0=hRhs*0;
            sol=solveKApeSymmetric(P,hL,f,hRhs,x0,y0,CtrlVar);
            varargout{k}=sol;

        else
            varargout{k}=P\varargin{k};
        end

        % If Cholesky available, I could do:
        % b=varargin{1} ;      [R,flag,P]=chol(MUA.M); xTest=P*(R\(R'\(P'*b))); x=MUA.M\b; [xTest(1:10) x(1:10)]

        % and if using vector:
        %  b=varargin{1} ;      [R,flag,p]=chol(MUA.M,"vector"); xTest2(p)=(R\(R'\(b(p)))); x=MUA.M\b; [xTest2(1:10) x(1:10)]



        if CtrlVar.Inverse.InfoLevel>=1000
            FindOrCreateFigure('I gradient') ;
            PlotMeshScalarVariable(CtrlVar,MUA,varargin{k}) ;
            hold on
            PlotMuaMesh(CtrlVar,MUA,nan,'w');
            title('Derivative Mesh Dependent')

            FindOrCreateFigure('P gradient') ;
            PlotMeshScalarVariable(CtrlVar,MUA,varargout{k}) ;
            hold on
            PlotMuaMesh(CtrlVar,MUA,nan,'w');
            title('Derivative Mesh Independent')
        end
    end
end


end
