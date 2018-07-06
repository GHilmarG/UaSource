function [gradJ1,RunInfo]=NewConjugatedGrad(dJdescent,dJdescentlast,gradJ0,CtrlVar,RunInfo)
    
    persistent iCount
    
    % [gradJ1,ConjGradAngle,teta]=NewConjugatedGrad(dJdescent,dJdescentlast,gradJ0,CtrlVar)
    % dJdescent and dJdecentlast are the current and the previous steepest descent directions
    % gradJ0 is previous gradient used in the minimisation and gradJ1 is the new gradient
    
    if isempty(iCount) ; iCount=0 ; end
    
    iCount=iCount+1;
    RunInfo.Inverse.ConjGradUpdate=RunInfo.Inverse.ConjGradUpdate+1;
    
    if isempty(CtrlVar)
        CtrlVar.ConjugatedGradientsRestartThreshold=0.5 ;
    elseif ~isfield(CtrlVar,'ConjugatedGradientsRestartThreshold')
        CtrlVar.ConjugatedGradientsRestartThreshold=0.5 ;
    end
    
    if isempty(CtrlVar)
        CtrlVar.ConjugatedGradientsUpdate='FR';
    elseif ~isfield(CtrlVar,'ConjugatedGradientsUpdate')
        CtrlVar.ConjugatedGradientsUpdate='FR';
    end
    
    
    d0=-dJdescentlast;
    d1=-dJdescent ;
    
    
    dd=d1'*d0/(norm(d1)*norm(d0));
    ddAngle=acosd(dd); % subsequent steepest descent directions should be close to 90 degrees angle to each other
    % dd=d1'*d0/(d0'*d0)
    
    if norm(dJdescent)<eps
        % return neg steepest descent direction
        ConjGradAngle=0 ; teta=0; gradJ1=d1;
        RunInfo.Inverse.ConjGradUpdate=0;
    elseif abs(ddAngle) < CtrlVar.ConjugatedGradientsRestartThreshold
        ConjGradAngle=0 ; teta=0; gradJ1=d1;
        if CtrlVar.Inverse.InfoLevel>=100
            fprintf(' resetting conjugated gradients because abs(ddAngle)=%g < CtrlVar.ConjugatedGradientsRestartThreshold=%g \n ',...
                abs(ddAngle),CtrlVar.ConjugatedGradientsRestartThreshold)
        end
        iCount=0;
        RunInfo.Inverse.ConjGradUpdate=0;
        
    else
        % seems that Polak-Ribi\ere is the best update
        switch upper(CtrlVar.ConjugatedGradientsUpdate)
            case 'FR'
                teta=d1'*d1/(d0'*d0);      % Fletcher-Reeves
            case 'PR'
                teta=d1'*(d1-d0)/(d0'*d0); % Polak-Ribi\`ere
            case 'HS'
                teta=d1'*(d1-d0)/(gradJ0'*(d1-d0)); % Hestenes-Stiefel
            case 'DY'
                teta=d1'*d1 /(gradJ0'*(d1-d0));   % Dai-Yan
            otherwise
                error('case not reckognized')
        end
        
        if teta<0
            
            if CtrlVar.Inverse.InfoLevel>=2
                fprintf(' resetting conjugated gradients because teta=%g<0. \n ',teta)
            end
            iCount=0;
            RunInfo.Inverse.ConjGradUpdate=0;
            teta=0;
            
        end
        
        if abs(ddAngle)<CtrlVar.ConjugatedGradientsRestartThreshold
            
            if CtrlVar.Inverse.InfoLevel>=2
                fprintf(' resetting conjugated gradients because angle between current and previous steepest descent direction')
                fprintf(' too far from 90 degrees. \n')
                fprintf(' abs(ddAngle)<%g. \n ',CtrlVar.ConjugatedGradientsRestartThreshold)
            end
            iCount=0;
            RunInfo.Inverse.ConjGradUpdate=0;
            teta=0;
            
        end
        
        
        
        
        gradJ1=d1 + teta * gradJ0 ; % new search direction
        ConjGradAngle=acosd(d1'*gradJ1/(norm(d1)*norm(gradJ1)));  % angle between current steepest descent and cc search direction
        
    end
    
    
    if CtrlVar.Inverse.InfoLevel>=2
        fprintf(' Conj. Grad update %s # %-i \n ',CtrlVar.ConjugatedGradientsUpdate,iCount)    
        fprintf(' teta=%g \n',teta)
        fprintf('angle between current and previous steepest descent directions is %-20.10g degrees \n',ddAngle)
        fprintf('     angle between conj. grad. and steepest descent directions is %-20.10g degrees \n',ConjGradAngle)
    end
    
end
