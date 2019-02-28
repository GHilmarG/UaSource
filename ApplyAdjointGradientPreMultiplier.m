
function varargout=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,varargin)



varargout=varargin;


switch CtrlVar.Inverse.AdjointGradientPreMultiplier
    
    case 'M'
        
        % The mass matrix has the dimensions area.
        Area=TriAreaTotalFE(MUA.coordinates,MUA.connectivity);
        
        if ~isfield(MUA,'M')
            MUA.M=MassMatrix2D1dof(MUA);
        end
        
        
        for k=1:numel(varargin)
            
            if ~isempty(varargin{k})
                varargout{k}=Area*(MUA.M\varargin{k});
                
                
                if CtrlVar.Inverse.InfoLevel>=1000
                    figure ; PlotMeshScalarVariable(CtrlVar,MUA,varargin{k}) ; title('Derivative Mesh Dependent')
                    figure ; PlotMeshScalarVariable(CtrlVar,MUA,varargout{k}) ; title('Derivative Mesh Independent')
                end
            end
        end
        
end
end
