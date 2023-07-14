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

    % R=Tint-Fext;
    % Tint=[Tx ; Ty ; Th] ;
    % Rint=[Fx ; Fy ; Fh] ;
    %
    % 
    %
    % Th=-dhdt 
    % Fh= a-dq/dx
    %
    %


    if ~isempty(L)
        
        frhs=-R-L'*luvh;   % Not sure why I put this minus here, but with the minus it becomes the right-hand side
        grhs=cuvh-L*[F1.ub;F1.vb;F1.h];
        
    else
        frhs=-R;
        grhs=[];
    end
    
    
%    d=[dub;dvb;dh]  ; % Newton step
%    D2=frhs'*d  ;

% The system I solve is K deltax = -R 
% and K=dR/dx
%
% D^2 = - R d = rhs d
%
%
    % Ru=-frhs(1:MUA.Nnodes) ;
    % Rv=-frhs(MUA.Nnodes+1:2*MUA.Nnodes)  ;
    % Rh=-frhs(2*MUA.Nnodes+1:3*MUA.Nnodes) ;
    % Rl=-grhs ; 
    % D2u=-Ru'*dub;
    % D2v=-Rv'*dvb;
    % D2h=-Rh'*dh;
    % D2l=-Rl'*dl; 
    % D2=full(D2u+D2v+D2h+D2l) ;  % For a feasable point, D2 must be positive for gamma=0, for the Newton direction to be a direction of decent for rWork
    
    D2=[frhs;grhs]'*[dub;dvb;dh;dl];
    rWork=full(D2^2);


    
    
    % rForce=ResidualCostFunction(CtrlVar,MUA,L,frhs,grhs,fext0,"-uvh-");
    % rForce=(frhs'*frhs+grhs'*grhs)/(fext0'*fext0+1000*eps);
    rForce=full([frhs;grhs]'*[frhs;grhs]./(fext0'*fext0+1000*eps));
    
    %% Testing TestIng
    % rForce=(R'*R)/(fext0'*fext0); 
    % rForce=full(rForce) ;
    %%

    switch CtrlVar.uvhMinimisationQuantity
        case "Force Residuals"
            r=rForce;
        case "Work Residuals"
            r=rWork;
        otherwise
            error("CalcCostFunctionNRuvh:CaseNotFound","case not found")
    end
    
    %     if nargout>=7
    %         Inodes=[];
    %         [ruv,rh]=CalcIncrementsNorm(CtrlVar,MUA,L,Inodes,F1.ub,dub,F1.vb,dvb,F1.h,dh);
    %     end
    
end

