
function varargout=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,varargin)



varargout=varargin;


switch CtrlVar.Inverse.AdjointGradientPreMultiplier
    
    case 'M'
        
        % The mass matrix has the dimensions area.

        H=MUA.M ;         
        
        % TestIng
        H=MUA.M+1e6*(MUA.Dxx+MUA.Dyy) ; 
        %
        
        
        for k=1:numel(varargin)
            
            if ~isempty(varargin{k})
                
                varargout{k}=MUA.Area*(H\varargin{k});
                
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
                    
                    FindOrCreateFigure('M gradient') ; 
                    PlotMeshScalarVariable(CtrlVar,MUA,varargout{k}) ; 
                    hold on
                    PlotMuaMesh(CtrlVar,MUA,[],'w');
                    title('Derivative Mesh Independent')
                end
            end
        end
        
end
end
