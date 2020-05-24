function [r,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,fext0)
    
    
    narginchk(15,15)
    nargoutchk(1,6)
    
    
    F1.ub=F1.ub+gamma*dub;
    F1.vb=F1.vb+gamma*dvb;
    F1.h=F1.h+gamma*dh;
    luvh=luvh+gamma*dl;
    
    
    CtrlVar.uvhMatrixAssembly.ZeroFields=false;
    CtrlVar.uvhMatrixAssembly.Ronly=true;
    
    
    [UserVar,RunInfo,R,~]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1);
    
    
    if ~isempty(L)
        
        Ruvh=-R-L'*luvh;   % Not sure why I put this minus here, but with the minus it becomes the right-hand side
        Rl=cuvh-L*[F1.ub;F1.vb;F1.h];
        
    else
        Ruvh=-R;
        Rl=[];
    end
    
    
%    d=[dub;dvb;dh]  ; % Newton step
%    D2=frhs'*d  ;

    Ru=-Ruvh(1:MUA.Nnodes) ;
    Rv=-Ruvh(MUA.Nnodes+1:2*MUA.Nnodes)  ;
    Rh=-Ruvh(2*MUA.Nnodes+1:3*MUA.Nnodes) ;
    
    % I= (F1.h <= (CtrlVar.ThickMin+1000*eps)); Ru(I)=0 ; Rv(I)=0 ; Rh(I)=0 ; 
 
    D2u=-Ru'*dub;
    D2v=-Rv'*dvb;
    D2h=-Rh'*dh;
    D2l=-Rl'*dl; 
    
    D2=D2u+D2v+D2h+D2l ;  % For a feasable point, D2 must be positive for gamma=0, for the Newton direction to be a direction of decent for rWork
    % Newton decrement: R' d
    
    rWork=D2^2; 

%    d=[dub;dvb;dh;dl]  ; D2=[Ruvh;Rl]'*d  ; rWork=D2^2 ;
    
    
    
    
    rForce=ResidualCostFunction(CtrlVar,MUA,L,Ruvh,Rl,fext0,"-uvh-");
    
    
    
    switch CtrlVar.uvhMinimisationQuantity
        case "Force Residuals"
            r=rForce;
        case "Work Residuals"
            r=rWork;
    end
    
%     if nargout>=7
%         Inodes=[]; 
%         [ruv,rh]=CalcIncrementsNorm(CtrlVar,MUA,L,Inodes,F1.ub,dub,F1.vb,dvb,F1.h,dh);
%     end
    
end

