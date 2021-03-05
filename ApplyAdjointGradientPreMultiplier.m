
function varargout=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,P,varargin)

narginchk(4,inf)


varargout=varargin;

if CtrlVar.Inverse.AdjointGradientPreMultiplier=="M"
    P=MUA.M/MUA.Area ;
elseif CtrlVar.Inverse.AdjointGradientPreMultiplier=="D"
    P=MUA.Dxx+MUA.Dyy  ;
elseif CtrlVar.Inverse.AdjointGradientPreMultiplier=="I" ...
        || CtrlVar.Inverse.AdjointGradientPreMultiplier=="Hanalytical" ...
        || isempty(P)
    varargout=varargin;
    return
end



for k=1:numel(varargin)
    
    if ~isempty(varargin{k})
        
        varargout{k}=P\varargin{k};
        
        %IgradientNorm=norm(varargin{k});
        %Mgradient=Area*(MUA.M\varargin{k});
        %MgradientNorm=norm(Mgradient) ;
        %varargout{k}=Mgradient*IgradientNorm/MgradientNorm;
        
        
        if CtrlVar.Inverse.InfoLevel>=1000
            FindOrCreateFigure('I gradient') ;
            PlotMeshScalarVariable(CtrlVar,MUA,varargin{k}) ;
            hold on
            PlotMuaMesh(CtrlVar,MUA,[],'w');
            title('Derivative Mesh Dependent')
            
            FindOrCreateFigure('P gradient') ;
            PlotMeshScalarVariable(CtrlVar,MUA,varargout{k}) ;
            hold on
            PlotMuaMesh(CtrlVar,MUA,[],'w');
            title('Derivative Mesh Independent')
        end
    end
end


end
