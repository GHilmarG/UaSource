


function [gamma,fgamma,LineSearchStagnated]=LineSearchNR(fdata,a,b,c,fa,fb,fc,f0,solrange,solbracket)
    
    InfoLevel=100;
    
    Slope=-f0;
    
    % f0 value at gamma=0
    % a<b<c on input, and fa fb fc function values at these points 
    
    LineSearchStagnated=0;
    fvector=[fa;fb;fc]; distancevector=[a;b;c];
    
    
    if fa+(fc-fa)*b/c < fb  && fa < fb
        fprintf(' curvature not positive and to the right of min , backtrack! \n ')  ;
        gamma=(a+b)/2; solbracket(2)=b;
    else
        gamma = parabolamin(a,b,c,fa,fb,fc); gamma=min([max([solrange(1) gamma]) solrange(2)]);
    end
    
      
    fgamma=fdata.func(gamma,fdata.arguments{:});
       
    
    fvector=[fvector ; fgamma ] ; distancevector=[distancevector ; gamma];
    [distancevector,isort]=sort(distancevector); fvector=fvector(isort) ;
    
    [smallestf,imin]=min(fvector); bestdistance=distancevector(imin);
    
    
    % gamma is the step length based on quadradic approximation
    % If the Newton-Raphson assumptions are fullfilled, we have gamma=1
    % usually this is the case exepct sometimes in the beginning and towards the end if residuals become comparable to
    % machine precision
    
    
    gammalast=-1000; tolLineSearch=0.001;
    
    ItNl=0;
    while smallestf > f0/2 || ItNl < 5  || (smallestf > f0 && ItNl <10 ) % only go in hear if f is not reduced by 50% in first atempt
        ItNl=ItNl+1;
        
        if imin>1 && imin<length(fvector)
            % fprintf(' case abc \n')
            solbracket=[distancevector(imin-1) ; distancevector(imin+1)] ;
            a=distancevector(imin-1) ; fa=fvector(imin-1);
            b=distancevector(imin)   ; fb=fvector(imin);
            c=distancevector(imin+1) ; fc=fvector(imin+1);
            gamma = parabolamin(a,b,c,fa,fb,fc); gamma=min([max([solrange(1) gamma]) solrange(2)]);
            % check if new gamma to close to old value
            if abs(gamma-gammalast)< (distancevector(2)-distancevector(1))/10;
                %fprintf(' bracketing ')
                if fa > fc
                    gamma=(a+b)/2;
                else
                    gamma=(b+c)/2;
                end
            end
        elseif imin==1
            fprintf(' imin=1 \n')
            dd=unique(sort(distancevector));
            if dd(2)< solbracket(2)
                solbracket(2)=dd(2) ;
            end
            gamma=(dd(1)+dd(2))/2;
        elseif imin==length(fvector)
            %fprintf(' imin=length(fvector) \n')
            if distancevector(imin) < solbracket(2) && distancevector(imin) > solbracket(1) ;
                solbracket(1)=distancevector(imin) ;
            elseif distancevector(imin-1) < solbracket(2) && distancevector(imin-1) > solbracket(1) ;
                solbracket(1)=distancevector(imin-1) ;
            end
            gamma=mean(solbracket);
        end
        
        if gamma < solbracket(1) || gamma > solbracket(2) ; fprintf(' gamma put within solbracet \n') ; gamma=mean(solbracket); end
        
        gamma=min([max([solrange(1) gamma]) solrange(2)]);

        fgamma=fdata.func(gamma,fdata.arguments{:});
        
        fprintf(' %g < %g < %g with fgamma % g , smallest f %g and f0 %g \n ',solbracket(1),gamma,solbracket(2),fgamma,smallestf,f0)
        
        fvector=[fvector ; fgamma ] ; distancevector=[distancevector ; gamma];
        [distancevector,isort]=sort(distancevector); fvector=fvector(isort) ;
        
        
        if solbracket(2)-solbracket(1) < tolLineSearch
            
            fprintf(' tolerance in bracketing reached, breaking \n')
            
            
            break
        end
        if abs(gamma-gammalast)< tolLineSearch && smallestf < f0
            fprintf(' change in gamma too small , breaking \n')
            break
        end
        
        if ItNl > 2 && gamma~=0 && smallestf < f0
            smallestvalues=unique(sort(fvector));
            if (smallestvalues(2)-smallestvalues(1))/f0 <0.01
                fprintf(' further reduction in function too small, breaking after %g iterations \n',ItNl)
                break
            end
        end
        
        [smallestf,imin]=min(fvector); bestdistance=distancevector(imin);
        gammalast=gamma;
        
        
    end
    
    
    
    if fgamma>smallestf && bestdistance>=solbracket(1) && bestdistance <= solbracket(2)
        gamma=bestdistance ; fgamma=smallestf ;
    end
    
    if gamma==0 ;
        fprintf(' gamma returned is zero ! Algorithim has stagnated \n ') ;
        LineSearchStagnated=1;
    end
    
  
    fgamma=fdata.func(gamma,fdata.arguments{:});
    
    figure(999) ; hold off ; plot(distancevector,fvector,'o') ; hold on ; plot(gamma,fgamma,'+r') ; plot(solbracket,[0;0],'xr') ;
    plot([0 ; 0.1],[f0;f0+Slope*0.1],'r')
    
    %        % step size selected by requiring residual force at the end of iteration to be orthogonal to (\Delta u, \Delta v)
    % 		R0=[du ; dv]'*R; R1=[du ; dv]'*R1;
    % 		a=R0/R1; asqr=sqrt((a/2)^2-a); ap=a/2+asqr;  am=a/2-asqr;
    % 		if a<0 ; gamma=ap ; else gamma=a/2; end
    % 		gamma=max([0.25 min([gamma 1])]);  % just in case quadradic appoximation is totally off
    % 		fprintf(' gamma %g R0 %g R1 %g ap %g am %g \n',gamma,R0,R1,ap,am)
    
    if InfoLevel > 0
        fprintf(' Line search gives gamma % g and Force residual for [u;v]+gamma [du;dv] %g , compared to %g for zero step size \n ',...
            gamma,fgamma, f0)
        
    end
end

