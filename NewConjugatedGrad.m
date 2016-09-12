function [gradJ1,ConjGradAngle,teta]=NewConjugatedGrad(dJdescent,dJdescentlast,gradJ0,CtrlVar)
    
    persistent iCount
    
    % [gradJ1,ConjGradAngle,teta]=NewConjugatedGrad(dJdescent,dJdescentlast,gradJ0,CtrlVar)
    % dJdescent and dJdecentlast are the current and the previous steepest descent directions
    % gradJ0 is previous gradient used in the minimisation and gradJ1 is the new gradient
    
    if isempty(iCount) ; iCount=0 ; end
    
    iCount=iCount+1;
    
    dd=dJdescent'*dJdescentlast/(norm(dJdescentlast)*norm(dJdescentlast));
    ddAngle=acosd(dd); % subsequent steepest descent directions should be close to 90 degrees angle to each other
    
    d0=-dJdescentlast;
    d1=-dJdescent ;
    
    
    % dd=d1'*d0/(d0'*d0)
    
    if norm(dJdescent)<eps
        % return neg steepest descent direction
        ConjGradAngle=0 ; teta=0; gradJ1=d1;
    elseif abs(dd) > CtrlVar.ConjugatedGradientsRestartThreshold
        ConjGradAngle=0 ; teta=0; gradJ1=d1;
        fprintf(' resetting conjugated gradients \n ')
        iCount=0;
    else
         % seems that Polak-Ribi\ere is the best update
        switch upper(CtrlVar.ConjugatedGradientsUpdate)
            case 'FR'
                teta=d1'*d1/(d0'*d0);      % Fletcher-Reeves
            case 'PR'
                teta=d1'*(d1-d0)/(d0'*d0); % Polak-Ribi\`ere
            case 'HS'
                teta=-d1'*(d1-d0)/(gradJ0'*(d1-d0)); % Hestenes-Stiefel
            case 'DY'
                teta=-d1'*d1 /(gradJ0'*(d1-d0));   % Dai-Yan
            otherwise
                error('case not reckognized')
        end

        gradJ1=d1 + teta * gradJ0 ; % new search direction
        ConjGradAngle=acosd(d1'*gradJ1/(norm(d1)*norm(gradJ1)));  % angle between current steepest descent and cc search direction
    end
    
    fprintf(' Conj. Grad update # %-i \n ',iCount)
    fprintf('angle between current and previous steepest descent directions is %-20.10g degrees \n',ddAngle)
    fprintf('     angle between conj. grad. and steepest descent directions is %-20.10g degrees \n',ConjGradAngle)
    
end
