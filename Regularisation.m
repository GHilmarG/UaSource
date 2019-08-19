function [R,dRdp,ddRddp,RegOuts]=Regularisation(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo)


RegOuts=[];
%%

% Add up field by field


[Areas,xEle,yEle,Area]=TriAreaFE(MUA.coordinates,MUA.connectivity);


%%
% I start by defining dpX, gsX and gaX, where X is either C or A and
%  dpX=X-X_{Prior}
%  gsX and gaX are the slope and amplitude regularization pre-factors
%

% C
if contains(lower(CtrlVar.Inverse.InvertFor),'c')
    
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


% b
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



if ~(CtrlVar.AGlenisElementBased   &&  CtrlVar.CisElementBased)
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
        %ddRdd=inv(Priors.CovC)/2/N;
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
        ddRCddpB=[];
    else
        RB=0;
        dRdB=[];
        
    end
    
    
    
    
    R=RAGlen+Rb+RC;
    dRdp=[dRdAGlen;dRdb;dRdC];
    
    
    
else  % Tikhonov regularization
    
    
    
    if isA
        if CtrlVar.AGlenisElementBased
            
            MA=gaA.^2*Areas/Area;
            npA=numel(dpA);
            NA=sparse(1:npA,1:npA,MA,npA,npA);
            
        else
            
            NA=(gsA.^2.*(Dxx+Dyy)+gaA.^2.*M)/Area;
            
        end
        RAGlen=dpA'*NA*dpA/2;
        dRdAGlen=(NA*dpA).*dAfactor;
    else
        RAGlen=0;
        dRdAGlen=[];
    end
    
    if isC
        if CtrlVar.CisElementBased
            
            MC=gaC.^2*Areas/Area;
            npC=numel(dpC);
            NC=sparse(1:npC,1:npC,MC,npC,npC);
            
        else
            
            NC=(gsC.^2.*(Dxx+Dyy)+gaC.^2.*M)/Area;
            
        end
        
        RC=dpC'*NC*dpC/2;
        dRdC=(NC*dpC).*dCfactor;
        
    else
        RC=0;
        dRdC=[];
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
        RB=dpB'*NB*dpB/2;
        dRdB=(NB*dpB).*dBfactor;
        
    else
        RB=0;
        dRdB=[];
    end
    
    [dRdAGlen,dRdb,dRdB,dRdC]=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,dRdAGlen,dRdb,dRdB,dRdC); 
        
    R=RAGlen+Rb+RB+RC;
    dRdp=[dRdAGlen;dRdb;dRdB;dRdC];
    ddRddp=[];
    
    
    
end


R=CtrlVar.Inverse.Regularize.Multiplier*R;
dRdp=CtrlVar.Inverse.Regularize.Multiplier*dRdp;
ddRddp=CtrlVar.Inverse.Regularize.Multiplier*ddRddp;


RegOuts.R=R;
RegOuts.dRdp=dRdp;
RegOuts.ddRddp=ddRddp;


RegOuts.RAGlen=RAGlen;
RegOuts.dRdAGlen=dRdAGlen;

RegOuts.RC=RC;
RegOuts.dRdC=dRdC;

RegOuts.Rb=Rb;
RegOuts.dRdb=dRdb;

RegOuts.RB=RB;
RegOuts.dRdB=dRdB;


