function [r,gamma,infovector,iarm,BacktrackInfo]=FindBestGammaBacktrackSSHEET(CtrlVar,r0,r1,F0,MUA,...
    AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dh,dt,Lh,lambdah,dlambdah)
         
    
    if CtrlVar.InfoLevelNonLinIt>=1000
        fprintf(CtrlVar.fidlog,...
            'FindBestGamma2DuvhBacktrack: on input r0=%-g \t r1=%-g \t r1/r0=%-g\n ',r0,r1,r1/r0) ;
    end
    
    BacktrackInfo.converged=1;  % true if it converged
    BacktrackInfo.gammaNaN=0;
    Slope0=-2*r0 ;  % using the inner product def
    gamma=1; r=r1; 
    gammab=1; rb=r1 ;
    gammac=1 ; rc=r1;
    
    beta=CtrlVar.NewtonBacktrackingBeta;  % Only accept step if reduction is larger than beta*StepSize
    target=max([CtrlVar.NewtonAcceptRatio*r0 CtrlVar.NLtol]);
    GammaMin=CtrlVar.BacktrackingGammaMin;
    
    
    iarm=0;
    iarmmax=CtrlVar.iarmmax;
    
    infovector=zeros(iarmmax+2,2)+NaN;
    infovector(1,1)=0 ; infovector(1,2)=r0;
    infovector(2,1)=1 ; infovector(2,2)=r1;
    
    cStatus=NaN ; pStatus=NaN;
    
    I=3;
    
    % when backtracking: gamma < gammab   and  gammac is the previous gammab with gammac> gammab
    % when extrapolating:  gammac=gamma ; gammab is the previous extrapolation value
    %                      gammab < gammac=gamma , rb> rc=r
    % initially r=r1 , rb=r1 and gamma=1, gammab=1
    
    
    %% possible initial extrapolation step
    
    if r> target && r1 < r0  && CtrlVar.LineSeachAllowedToUseExtrapolation
        ExtrapolationStep=true ;
        if CtrlVar.InfoLevelNonLinIt>=2
            fprintf(CtrlVar.fidlog,' Extrapolation flag set to true in line-search uv \n ');
        end
    else
        ExtrapolationStep=false;
    end
    
    while ExtrapolationStep && iarm<=10
        % In an extrapolation step I do not use the info on slope at 0
        if iarm==0
            gammaTest=1.5;  % I need a reasonable initial guess
        elseif iarm==1
            [gammaTest,cstatus] = parabolamin(0,1,gammac,r0,r1,rc);
        elseif iarm==2
            [gammaTest,cstatus] = parabolamin(1,gammab,gammac,r1,rb,rc);
        else
            [gammaTest,cstatus] = parabolamin(gammaa,gammab,gammac,ra,rb,rc);
        end
        
        if gammaTest > 1.1*gammac  %  suggests extrapolation and new gamma not too close to previous one
            gamma=gammaTest;
            if gamma>2* gammac ; gamma=2*gammac ; end  % guard against wild extrapolation
            
            r=CalcCostFunctionSSHEET(CtrlVar,gamma,dh,MUA,AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dt,Lh,lambdah,dlambdah,F0);
            %[r,ruv,rh]=CalcCostFunctionNRuvh(gamma,du,dv,dh,u,v,h,S,B,b,u0,v0,h0,as0,ab0,as1,ab1,dt,AGlen,n,C,m,coordinates,connectivity,nip,alpha,rho,rhow,g,MeshProp,Boundary,CtrlVar,F0,L,lambda,dlambda);
            infovector(I,1)=gamma ; infovector(I,2)=r; I=I+1;
            
            if isnan(r) ;
                fprintf(CtrlVar.fidlog,' In line search (extrapolation phase), cost function is nan \n ');
                error(' nan ')
            end
            
            if r < rc
                if r>0.99*rc ;  ExtrapolationStep=false; end  % don't bother with extrapolation if reduction is less than 1%
                iarm=iarm+1;
                gammaa=gammab ; ra=rb ;
                gammab=gammac ; rb=rc ;
                gammac=gamma  ; rc=r;
                target=max([ r0-beta*Slope0*gamma, CtrlVar.NLtol])  ; % Armijo criteria
                if CtrlVar.InfoLevelNonLinIt>=2
                    fprintf(CtrlVar.fidlog,'E: iarm=%-i \t gammab=%-g \t gammac=%-g \t r0=%-g \t r1=%-g \t r=%-g \t  rc=%-g \t target=%-g \t r/rtarget=%-g \n',iarm,gammab,gammac,r0,r1,r,rc,target,r/target);
                end
            else
                ExtrapolationStep=false;
            end
        else
            ExtrapolationStep=false;
        end
    end
    
    
    
    
    %% backtracking step
    while r >  target && iarm<=iarmmax && gamma > GammaMin 
        
        iarm=iarm+1;
        if iarm==1  % parabolic backtracking step
            gamma=-Slope0/2/(rb-r0-Slope0);
            if gamma > 0.8*gammab ; gamma=0.8*gammab ; elseif gamma < 0.4*gammab ; gamma=0.4*gammab;end
        else % cubic backtracking step
            [gamma,cStatus]=CubicFit(Slope0,r0,rb,rc,gammab,gammac);
            if cStatus==1 ;
                fprintf(CtrlVar.fidlog,'Cubic Fit returns status 1 with gamma=%-g \n ',gamma);
                [gamma,pStatus] = parabolamin(0,gammab,1,r0,rb,r1);
                if pStatus==1
                    fprintf(CtrlVar.fidlog,'parabolamin returns status 1 with gamma%-g \n ',gamma);
                end
            end
            if gamma > 0.8*gammab ; gamma=0.8*gammab ; elseif gamma < 0.2*gammab ; gamma=0.2*gammab;end
        end
        
        if isnan(gamma) ;
            warning('FindBestGammaBacktrackSSHEET:gammaNaN',' gamma is NaN ') ;
            fprintf(CtrlVar.fidlog,' gamma in backtracking is NaN \n ');
            save TestSave Slope0 r0 r1 rb rc gammab gammac cStatus pStatus
            BacktrackInfo.gammaNaN=1;
            break
        end
        
        r=CalcCostFunctionSSHEET(CtrlVar,gamma,dh,MUA,AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dt,Lh,lambdah,dlambdah,F0); 
        %[r,ruv,rh]=CalcCostFunctionNRuvh(gamma,du,dv,dh,u,v,h,S,B,b,u0,v0,h0,as0,ab0,as1,ab1,dt,AGlen,n,C,m,coordinates,connectivity,nip,alpha,rho,rhow,g,MeshProp,Boundary,CtrlVar,F0,L,lambda,dlambda);
        infovector(I,1)=gamma ; infovector(I,2)=r; I=I+1;
        
        if isnan(r) ;
            fprintf(CtrlVar.fidlog,' In line search (backtracking phase), cost function is nan \n ');
            error(' nan ')
        end
        
        
        if iarm==1  &&  rb< r ;
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
        
        
        target= max([r0-beta*Slope0*gamma,CtrlVar.NLtol])  ; % Armijo criteria
        if CtrlVar.InfoLevelNonLinIt>=2
            fprintf(CtrlVar.fidlog,'B: iarm=%-i \t gammab=%-g \t gammac=%-g \t r0=%-g \t r1=%-g \t r=%-g \t  rc=%-g \t target=%-g \t r/rtarget=%-g \n',iarm,gammab,gammac,r0,r1,r,rc,target,r/target);
        end
    end
    
    
    
    if iarm==0 && r> 0.1*r0 && r > CtrlVar.NLtol
        %% it is so cheap to calculate the cost function that I do so even if first step has been accepted (here I'm allowing for some extrapolation beyond full Newton step)
        
        gammaTest=-Slope0/2/(rb-r0-Slope0);
        if gammaTest > 1.2*gammab ; gammaTest=1.2*gammab ; elseif gammaTest < 0.1*gammab ; gammaTest=0.1*gammab;end
        rTest=CalcCostFunctionSSHEET(CtrlVar,gamma,dh,MUA,AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dt,Lh,lambdah,dlambdah,F0); 
        %[rTest,ruvTest,rhTest]=CalcCostFunctionNRuvh(gammaTest,du,dv,dh,u,v,h,S,B,b,u0,v0,h0,as0,ab0,as1,ab1,dt,AGlen,n,C,m,coordinates,connectivity,nip,alpha,rho,rhow,g,MeshProp,Boundary,CtrlVar,F0,L,lambda,dlambda);
        infovector(I,1)=gammaTest ; infovector(I,2)=rTest; I=I+1;
        if rTest<r ; r=rTest ;  gamma=gammaTest ; iarm=iarm+1; end
        
    elseif iarm==1 && r> 0.1*r0 && BacktrackInfo.gammaNaN==0 && r > CtrlVar.NLtol
        % another case that is sometimes worthwile investigating is if first backtracking step can be improved
        
        gammaTest=gamma/2;
        rTest=CalcCostFunctionSSHEET(CtrlVar,gamma,dh,MUA,AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dt,Lh,lambdah,dlambdah,F0); 
        %[rTest,ruvTest,rhTest]=CalcCostFunctionNRuvh(gammaTest,du,dv,dh,u,v,h,S,B,b,u0,v0,h0,as0,ab0,as1,ab1,dt,AGlen,n,C,m,coordinates,connectivity,nip,alpha,rho,rhow,g,MeshProp,Boundary,CtrlVar,F0,L,lambda,dlambda);
        infovector(I,1)=gammaTest ; infovector(I,2)=rTest; I=I+1;
        if rTest<r ; r=rTest ; gamma=gammaTest ;  iarm=iarm+1; end
        
        [~,ind]=sort(infovector(:,1)); infovector(:,1)=infovector(ind,1); infovector(:,2)=infovector(ind,2);
        
        % cubic fit trough the values at 0, 1, and gamma
        gammaTest=CubicFit(Slope0,r0,infovector(2,2),1,infovector(2,1),r1);
        if gammaTest > 0.9 ; gammaTest=0.9 ; elseif gammaTest < 0.01 ; gammaTest=0.01 ;end
        rTest=CalcCostFunctionSSHEET(CtrlVar,gamma,dh,MUA,AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dt,Lh,lambdah,dlambdah,F0); 
        infovector(I,1)=gammaTest ; infovector(I,2)=rTest; I=I+1;
        if rTest<r ; r=rTest ; gamma=gammaTest ;
            fprintf(CtrlVar.fidlog,' cubic fit improved \n ');
        end
        
    end     
    infovector=infovector(1:I-1,:);
    
    if iarm>iarmmax
        warning('FindBestGammaBacktrackSSHEET:TooManyBacktrackingSteps',' maximum number of backtracking steps reached ')
        fprintf(CtrlVar.fidlog,' Maximum number of backtracking steps reached! ') ;
    end
    
    if r > r0 && r > CtrlVar.NLtol;
        fprintf(CtrlVar.fidlog,' Backtracking step did not reduce residual! r0=%-g and r=%-g \n ',r0,r) ;
        BacktrackInfo.converged=0;
    end
    
end

