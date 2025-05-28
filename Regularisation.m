


function [R,dRdp,ddRdpp,RegOuts]=Regularisation(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo)

%%
% Calculates the regularization term R, and the gradient and the Hessian of R with respect to p.
%
%
%
% This is a fairly simple thing to do as the regularization term is an explicit function of p, and the Hessian calculation can
% be done exactly.
%
% However, there are quite a few cases to consider...
%
% For each variable (A,B,C) the regularization term typically has the form
%
%   R= ( ga.^2* Ra + gs.^2 * Rs ) /(2*Area)
%
% where ga and gs are regularization parameters
%
% and 
%
%  Ra = dp'*M*dp
%
%  Rs= dp'*(Dxx_Dyy)*dp
%
% where dp = p-pPrior, and M is the mass matrix and Dxx and Dyy the stiffness matrices.
%
%
% Note: Although this focus is on what is generally considered to be regularization terms, these terms can also be thought of as
% those where the inverted fields (e.g. B C A) explicitly enter the cost function.
%
% Here we, for example, also include the deviation of B from direct measurements. This could equally be referred to as a misfit term
% as well. This B 'misfit' term has the form:
%
% $$ \int (B-B_{\mathrm{meas}}) \, C^{-1} \, (B-B_{\mathrm{meas}}) \, dA $$
%
% By defining:
%
% $$ \tilde{B} = (B-B_{\mathrm{meas}})/B_{\epsilon} $$
% 
% This simply has the form
%
% $$ \tilde{B} \, M \tilde{B} $$
%
%%


if nargout > 3
    RegOuts=[];
    RegOuts.RAs=nan  ; RegOuts.RAa=nan;
    RegOuts.RCs=nan  ; RegOuts.RCa=nan;
   
end
%%

% Add up field by field

RAGlen=0;
dRdAGlen=[];
ddRdAA=[];

RC=0;
dRdC=[];
ddRdCC=[];


Area=MUA.Area; 


%%
% I start by defining dpX, gsX and gaX, where X is either C or A and
%  dpX=X-X_{Prior}
%  gsX and gaX are the slope and amplitude regularization pre-factors
%

% C
if contains(lower(CtrlVar.Inverse.InvertFor),'c')  % this includes both c and logc inversion
    
    isC=1;
    if contains(lower(CtrlVar.Inverse.Regularize.Field),'logc')
        
        % regularize log10(C)
        dpC=log10(F.C)-log10(Priors.C);
        pPriorCovC=Priors.CovC;
        gsC=CtrlVar.Inverse.Regularize.logC.gs;
        gaC=CtrlVar.Inverse.Regularize.logC.ga;
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
            dCfactor=1;
        else
            dCfactor=1./F.C/log(10); % d (dpC)/dC= 1/(log(10) C)
        end
        
    else
        % regularize C
        dpC=F.C-Priors.C;
        pPriorCovC=Priors.CovC;
        gsC=CtrlVar.Inverse.Regularize.C.gs;
        gaC=CtrlVar.Inverse.Regularize.C.ga;
        if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
            dCfactor=log(10)*F.C;   % gradient must be with respect to logC, but regularisation is on C
        else
            dCfactor=1;
        end
    end
else
    
    isC=0;
    dpC=0;
    dCfactor=0;
    pPriorCovC=1;
    gsC=0;
    gaC=0;
    
end

% AGlen
if contains(lower(CtrlVar.Inverse.InvertFor),'aglen')
    
    isA=1;
    if contains(lower(CtrlVar.Inverse.Regularize.Field),'logaglen')
        
        % regularize log10(AGlen)
        
        dpA=log10(F.AGlen)-log10(Priors.AGlen);
        pPriorCovA=Priors.CovAGlen;
        gsA=CtrlVar.Inverse.Regularize.logAGlen.gs;
        gaA=CtrlVar.Inverse.Regularize.logAGlen.ga;
        if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            dAfactor=1;
        else
            dAfactor=1./F.AGlen/log(10);   % gradient must be with respect to A, but regularisation is on logA
        end
        
    else % regularize A
        
        dpA=F.AGlen-Priors.AGlen;
        pPriorCovA=Priors.CovAGlen;
        gsA=CtrlVar.Inverse.Regularize.AGlen.gs;
        gaA=CtrlVar.Inverse.Regularize.AGlen.ga;
        if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            dAfactor=log(10)*F.AGlen;   % gradient must be with respect to logA, but regularisation is on A
        else
            dAfactor=1;
        end
    end
else
    
    isA=0;
    dpA=0;
    dAfactor=0;
    pPriorCovA=1;
    gsA=0;
    gaA=0;
    
end

% b
if contains(CtrlVar.Inverse.InvertFor,'-b-')
    error('fdsa')
    isb=1;
    dpb=F.b-Priors.b;
    pPriorCovb=Priors.Covb;
    gsb=CtrlVar.Inverse.Regularize.b.gs;
    gab=CtrlVar.Inverse.Regularize.b.ga;
    dbfactor=1;
    
else
    
    isb=0;
    dpb=0;
    dbfactor=0;
    pPriorCovb=1;
    gsb=0;
    gab=0;
    
end


% B
if contains(CtrlVar.Inverse.InvertFor,'-B-')
    
    isB=1;
    dpB=F.B-Priors.B;
    pPriorCovB=Priors.CovB;
    gsB=CtrlVar.Inverse.Regularize.B.gs;
    gaB=CtrlVar.Inverse.Regularize.B.ga;
    dBfactor=1;
    
else
    
    isB=0;
    dpB=0;
    dBfactor=0;
    pPriorCovB=1;
    gsB=0;
    gaB=0;
    
end

if ~isfield(MUA,'M')
    M=MassMatrix2D1dof(MUA);
else
    M=MUA.M;
end

if ~isfield(MUA,'Dxx')
    [Dxx,Dyy]=StiffnessMatrix2D1dof(MUA);
else
    Dxx=MUA.Dxx ;
    Dyy=MUA.Dyy;
end




%% Now dpX, gsX and gaX have all be defined

% Now defining R, dRdp, ddRddp 
if contains(lower(CtrlVar.Inverse.Regularize.Field),'cov')  % Bayesian regularization
    
    % R= (C-C_prior)' CC^{-1} (C-C_prior)  / (2N)
    
    if isA
        npA=numel(dpA);
        temp=pPriorCovA\dpA;
        RAGlen=dpA'*temp/(2*npA)   ;
        dRdAGlen=temp/npA;
        %ddRdAA=inv(Priors.CovC)/2/N;
        ddRAddpA=[];
    else
        RAGlen=0;
        dRdAGlen=[];
        
    end
    
    if isC
        npC=numel(dpC);
        temp=pPriorCovC\dpC;
        RC=dpC'*temp/(2*npC)   ;
        dRdC=temp/npC;
        %ddRdd=inv(Priors.CovC)/2/N;
        ddRCddpC=[];
    else
        RC=0;
        dRdC=[];
        
    end
    
    if isb
        npb=numel(dpb);
        temp=pPriorCovb\dpb;
        Rb=dpb'*temp/(2*npb)   ;
        dRdb=temp/npb;
        %ddRdd=inv(Priors.CovC)/2/N;
        ddRCddpb=[];
        error('fdsa')
    else
        Rb=0;
        dRdb=[];
        
    end
    
    
    if isB
        npB=numel(dpB);
        temp=pPriorCovB\dpB;
        RB=dpB'*temp/(2*npB)   ;
        dRdB=temp/npB;
        %ddRdd=inv(Priors.CovC)/2/N;
        ddRddpB=[];
    else
        RB=0;
        dRdB=[];
        
    end
    
    
    
    
    R=RAGlen+Rb+RC;
    dRdp=[dRdAGlen;dRdb;dRdC];
    
    
    
else  % Tikhonov regularization
    
    
    
    if isA
        
        

        NA=(gsA.^2.*(Dxx+Dyy)+gaA.^2.*M)/Area;
        %RAGlen=dpA'*NA*dpA/2;
        dRdAGlen=(NA*dpA).*dAfactor;
        
        
        RAs= dpA'*(Dxx+Dyy)*dpA   / (2*Area);
        RAa= dpA'    *M    *dpA   /(2*Area);
        RAGlen=gsA.^2*RAs+gaA.^2*RAa; 
     
        RegOuts.RAs=RAs  ; RegOuts.RAa=RAa;


        if contains(CtrlVar.Inverse.MinimisationMethod,"Hessian")
            
            if contains(CtrlVar.Inverse.Hessian,"RHA=E")
                ddRdAA=NA.*dAfactor;
            elseif contains(CtrlVar.Inverse.Hessian,"RHA=M")
                ddRdAA=MUA.M/MUA.Area;
            elseif contains(CtrlVar.Inverse.Hessian,"RHA=I") || contains(CtrlVar.Inverse.Hessian,"RHA=1")
                N=MUA.Nnodes;
                ddRdAA=speye(N,N);
            elseif contains(CtrlVar.Inverse.Hessian,"RHA=0") || contains(CtrlVar.Inverse.Hessian,"RHA=O")
                N=MUA.Nnodes;
                ddRdAA=sparse(N,N);
            else
                fprintf(" CtrlVar.Inverse.Hessian=%s, is incorrect.\n",CtrlVar.Inverse.Hessian)
                error("Regularisation:IncorrectInputs"," case not found ")
            end
        end
    end
    
    if isC
        
    
        
        % RCs should always be positive. However, I discovered that it can happen that the smallest eigenvalue is slightly
        % negative!!! This must be due to numerical rounding errors when assembling Dxx and Dyy. I for example found a case where the
        % two smallest eigenvalues of Dyy were -1.14405445408737e-16 and  -8.99887803162969e-17. One approach of dealing with this
        % would be to add eps to the diagonal of Dxx and Dyy.
     
        Ieps=sparse(1:MUA.Nnodes,1:MUA.Nnodes,eps);
        Dxx=Dxx+Ieps ; Dyy=Dyy+Ieps ; 

        NC=(gsC.^2.*(Dxx+Dyy)+gaC.^2.*M)/Area;
        %RC=dpC'*NC*dpC/2;
        dRdC=(NC*dpC).*dCfactor;

        RCs= dpC'*(Dxx+Dyy)*dpC   / (2*Area);
        RCa= dpC'    *M    *dpC   /(2*Area);
        RC=gsC.^2*RCs+gaC.^2*RCa; 
        

        RegOuts.RCs=RCs  ; RegOuts.RCa=RCa;

        if contains(CtrlVar.Inverse.MinimisationMethod,"Hessian")
            if contains(CtrlVar.Inverse.Hessian,"RHC=E")
                
                ddRdCC=NC.*dCfactor;
            elseif contains(CtrlVar.Inverse.Hessian,"RHC=M")
                ddRdCC=MUA.M/MUA.Area;
            elseif contains(CtrlVar.Inverse.Hessian,"RHC=I") || contains(CtrlVar.Inverse.Hessian,"RHC=1")
                N=MUA.Nnodes;
                ddRdCC=speye(N,N);
            elseif contains(CtrlVar.Inverse.Hessian,"RHC=O") || contains(CtrlVar.Inverse.Hessian,"RHC=0")
                N=MUA.Nnodes;
                ddRdCC=sparse(N,N);
            else
                
                fprintf(" CtrlVar.Inverse.Hessian=%s, is incorrect.\n",CtrlVar.Inverse.Hessian)
                error("Regularisation:IncorrectInputs"," case not found ")
                
            end
        end
        
    end
    
    
    if isb   %  b
        
        Nb=(gsb.^2.*(Dxx+Dyy)+gab.^2.*M)/Area;
        Rb=dpb'*Nb*dpb/2;
        dRdb=(Nb*dpb).*dbfactor;
        error('fdsa')
    else
        Rb=0;
        dRdb=[];
    end



    if isB   %  B

        NB=(gsB.^2.*(Dxx+Dyy)+gaB.^2.*M)/Area;
        RB=dpB'*NB*dpB/2;               %       R: Regularisation term for B (a scalar)
        dRdB=(NB*dpB).*dBfactor;        %   dR/dB:  (a vector)
        ddRdBB=NB.*dBfactor;            % exact, or simply the correct, Hessian of the regularization term
        % To do: I could add "RHB=E" to CtrlVar.Inverse.Hessian. Right now I do the exact (E) Hessian evaluation here.


        if ~isempty(Meas.B)  &&  ~isempty(Meas.BCov)  &&    isdiag(Meas.BCov)

            % Adding a cost term giving the deviation of inverted B from direct measurements of B. This has the same form as a data
            % misfit term used for velocities and dh/dt. But here this is applied to the inverted field.
            

            Berr=sqrt(spdiags(Meas.BCov));
            Bres=(F.B-Meas.B)./Berr;
            RBmeas=full(Bres'*MUA.M*Bres)/2/Area;
            dRdBmeas=(MUA.M*Bres)./Berr/Area;
            ddRdBmeasBmeas=(MUA.M)./Berr/Area;

            RB=RB+RBmeas;
            dRdB=dRdB+dRdBmeas;
            ddRdBB=ddRdBmeasBmeas;

        end
        %

    else
        RB=0;
        dRdB=[];
        ddRdBB=[];
    end

    
    % if CtrlVar.Inverse.MinimisationMethod contains "Hessian", then the pre-multipler is simply I, so this has no effect.
    dRdAGlen=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,ddRdAA,dRdAGlen);
    dRdC=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,ddRdCC,dRdC);
    dRdB=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,ddRdBB,dRdB);
    
    R=RAGlen+RB+RC;
    dRdp=[dRdAGlen;dRdB;dRdC];
    
    
    [Am,An] = size(ddRdAA);
    [Bm,Bn] = size(ddRdBB);
    [Cm,Cn] = size(ddRdCC);
    ddRdpp = spalloc(Am+Bm+Cm,An+Bn+Cn,nnz(ddRdAA)+nnz(ddRdCC)+nnz(ddRdBB));
    ddRdpp(1:Am,1:An) = ddRdAA;
    ddRdpp(Am+1:Am+Bm,An+1:An+Bn) = ddRdBB;
    ddRdpp(Am+Bm+1:Am+Bm+Cm,An+Bn+1:An+Bn+Cn) = ddRdCC;
    
    
end


R=CtrlVar.Inverse.Regularize.Multiplier*R;
dRdp=CtrlVar.Inverse.Regularize.Multiplier*dRdp;
ddRdpp=CtrlVar.Inverse.Regularize.Multiplier*ddRdpp;

if nargout > 3
    RegOuts.R=R;
    RegOuts.dRdp=dRdp;
    RegOuts.ddRddp=ddRdpp;


    RegOuts.RAGlen=RAGlen;
    RegOuts.dRdAGlen=dRdAGlen;

    RegOuts.RC=RC;
    RegOuts.dRdC=dRdC;

    RegOuts.Rb=Rb;
    RegOuts.dRdb=dRdb;

    RegOuts.RB=RB;
    RegOuts.dRdB=dRdB;
end

if R< 0
    fprintf("Regularisation.m : R is negative \n")
end


end