function [UserVar,r,gamma,infovector,BacktrackInfo] = FindBestGamma2Dbacktracking(UserVar,CtrlVar,MUA,F,F0,r0,r1,L,l,cuv,dub,dvb,dl)
      
%FindBestGamma2Dbacktracking(UserVar,CtrlVar,MUA,F0,r0,r1,s,S,B,h,ub,dub,vb,dvb,uo,vo,AGlen,n,C,m,alpha,rho,rhow,g,L,l,dl,cuv)
   

    if ~isfield(CtrlVar,'NLtol')
        CtrlVar.NLtol=1e-15; 
    end
    
    Slope0=-2*r0 ;  % using the inner product def
    gamma=1; r=r1;
    gammab=1; rb=r1 ; 
    gammac=1 ; rc=r1;
    
  
    
    beta=CtrlVar.NewtonBacktrackingBeta;  % Only accept step if reduction is larger than beta*StepSize
    target=CtrlVar.NewtonAcceptRatio*r0;  % the initial criteria is a given fractional reduction, in the backtracking step I then use the Amarijo criteria
    GammaMin=CtrlVar.BacktrackingGammaMin;
    
    
    
    iarm=0;
    iarmmax=CtrlVar.iarmmax;
    
    
    BacktrackInfo.Converged=1;
    BacktrackInfo.gammaNaN=0;
    BacktrackInfo.iarm=iarm;
    
    infovector=zeros(iarmmax+2,2)+NaN;
    infovector(1,1)=0 ; infovector(1,2)=r0;
    infovector(2,1)=1 ; infovector(2,2)=r1;
    
%     
%     
%     if r1<CtrlVar.NLtol
%         return
% %    elseif r0<CtrlVar.NLtol
% %        gamma=0;
% %        r=r0;
%     end
%     
%     
    cStatus=NaN ; pStatus=NaN;
    
    I=3;
    
    % when backtracking: gamma < gammab   and  gammac is the previous gammab with gammac> gammab
    % when extrapolating:  gammac=gamma ; gammab is the previous extrapolation value
    %                      gammab < gammac=gamma , rb> rc=r
    % initially r=r1 , rb=r1 and gamma=1, gammab=1 
    
    %% possible initial extrapolation step
    if r> target && r1 < r0  && CtrlVar.LineSearchAllowedToUseExtrapolation
        ExtrapolationStep=true ;
        if CtrlVar.InfoLevelNonLinIt>=2
            fprintf(CtrlVar.fidlog,' Extrapolation flag set to true in line-search uv \n ');
        end
    else
        ExtrapolationStep=false;
    end
        
    while ExtrapolationStep && iarm<=20
        % In an extrapolation step I do not use the info on slope at 0
        if iarm==0 
            gammaTest=1.5;  % I need a reasonable initial guess
        elseif iarm==1 
            gammaTest = parabolamin(0,1,gammac,r0,r1,rc);
        elseif iarm==2
            gammaTest = parabolamin(1,gammab,gammac,r1,rb,rc);
        else
             gammaTest = parabolamin(gammaa,gammab,gammac,ra,rb,rc);
        end
        
        if gammaTest > 1.1*gammac  %  suggests extrapolation and new gamma not too close to previous one
            gamma=gammaTest;
            if gamma>2* gammac ; gamma=2*gammac ; end  % guard against wild extrapolation
                
            
           
            r = CalcCostFunctionNR(UserVar,[],CtrlVar,MUA,gamma,F,F0,L,l,cuv,dub,dvb,dl); 
           
            infovector(I,1)=gamma ; infovector(I,2)=r; I=I+1;
            
            if isnan(r) 
                fprintf(CtrlVar.fidlog,' In line search (extrapolation phase), cost function is nan \n ');
                error(' nan ')
            end
            
            if r < rc  
                if r>0.99*rc ;  ExtrapolationStep=false; end  % don't bother with extrapolation if reduction is less than 1%
                iarm=iarm+1;
                gammaa=gammab ; ra=rb ;
                gammab=gammac ; rb=rc ;
                gammac=gamma  ; rc=r;
                target= r0-beta*Slope0*gamma  ; % Armijo criteria
                if CtrlVar.InfoLevelNonLinIt>=2
                    fprintf(CtrlVar.fidlog,'E: iarm=%-i \t gammab=%-g \t gammac=%-g \t r0=%-g \t r1=%-g \t r=%-g \t  rb=%-g \t target=%-g \t r/rtarget=%-g \n',iarm,gammab,gammac,r0,r1,r,rb,target,r/target);
                end
                
                
            else
                ExtrapolationStep=false;
            end
        else
            ExtrapolationStep=false;
        end
        
    end
        
    iarmExtrapolation=iarm;
    iarm=0;
    
    %% backtracking
    while r >  target && iarm<=iarmmax && gamma > GammaMin  % && r > CtrlVar.NLtol
    
        
        iarm=iarm+1;
        if iarm==1  % parabolic fit
            gamma=-Slope0/2/(rb-r0-Slope0);
            if gamma > 0.8*gammab ; gamma=0.8*gammab ; elseif gamma < 0.2*gammab ; gamma=0.2*gammab;end
        else % cubic fit
            gamma=CubicFit(Slope0,r0,rb,rc,gammab,gammac); 
            if gamma > 0.5*gammab ; gamma=0.5*gammab ; elseif gamma < 0.2*gammab ; gamma=0.2*gammab;end
        end
             
        if isnan(gamma) 
            warning('FindBestGamma2DuvhBacktrack:gammaNaN',' gamma is NaN ') ;
            fprintf(CtrlVar.fidlog,' gamma in backtracking is NaN \n ');
            save TestSave Slope0 r0 r1 rb rc gammab gammac cStatus pStatus
            BacktrackInfo.gammaNaN=1;
            break
        end
        
       
        r = CalcCostFunctionNR(UserVar,[],CtrlVar,MUA,gamma,F,F0,L,l,cuv,dub,dvb,dl); 
       
        infovector(I,1)=gamma ; infovector(I,2)=r; I=I+1;
   
        if isnan(r) 
            fprintf(CtrlVar.fidlog,' In line search (backtracking phase), cost function is nan \n ');
            error(' nan ')
        end
        
        if iarm==1  &&  rb< r 
            rc=rb; gammac=gammab;
            rTemp=r; gammaTemp=gamma;
            r=rb ; gamma=gammab ;        % Because I start with a different exit criterion, I must
            rb=rTemp ; gammab=gammaTemp ;   % check if initial estimate is actually better than first improvement
            if CtrlVar.InfoLevelNonLinIt>10
                fprintf(CtrlVar.fidlog,' Newton step better than first backtracking step \n');
            end
        else
            gammac=gammab; gammab=gamma; rc=rb; rb=r; % if backtracking then this is correct
        end
        
        target= r0-beta*Slope0*gamma  ; % Armijo criteria
        if CtrlVar.InfoLevelNonLinIt>1
            fprintf(CtrlVar.fidlog,'B: iarm=%-i \t gammab=%-g \t gammac=%-g \t r0=%-g \t rb=%-g \t  rc=%-g \t target=%-g \t r/rtarget=%-g \n',iarm,gammab,gammac,r0,rb,rc,target,r/target);
        end
    end
    
    
    

    %% it is so cheap to calculate the cost function that I do so even if first step has been accepted 
        
    % another case that is sometimes worthwile investigating is if first backtracking step was not agressive enough
    if iarm==1 && r> 0.1*r0 && r> CtrlVar.NLtol
        gammaTest=gamma/2;
        
        
        rTest = CalcCostFunctionNR(UserVar,[],CtrlVar,MUA,gammaTest,F,F0,L,l,cuv,dub,dvb,dl); 

        infovector(I,1)=gammaTest ; infovector(I,2)=rTest; I=I+1;
        if rTest<r ; r=rTest ; gamma=gammaTest ; end
    end
    
    
    %%  
    
        
    if iarm>iarmmax
        warning('MATLAB:FindBestGamma2Dbacktracking:TooManyBacktrackingSteps',' maximum number of backtracking steps reached ')
        
    end
    
    if r>r0 && r>CtrlVar.NLtol % if everyting else fails, allow some increase, possibly it must get out of a local minimum
        fprintf(CtrlVar.fidlog,'FindBestGamma: residual increased from %g to %g, but still returning the new value \n !',r0,r);
        BacktrackInfo.Converged=0;    
    end
    
    infovector=infovector(1:I-1,:);
    
    if isnan(gamma) ; save TestSave ; error(' gamma is nan \n') ; end
    
    BacktrackInfo.iarm=iarm+iarmExtrapolation;    
    
end

