function F=p2F(CtrlVar,MUA,p,F,Meas,Priors)

persistent GLgeo GLnodes GLele



narginchk(6,6)

NA=numel(F.AGlen);
Nb=numel(F.b);
NB=numel(F.B);
NC=numel(F.C);

IA1=0 ; IA2=0;
Ib1=0 ; Ib2=0;
IB1=0 ; IB2=0;
IC1=0 ; IC2=0;
Itotal=0;

isA=false ; isb=false ; isB=false ; isC=false;

if contains(CtrlVar.Inverse.InvertForField,'A')
    IA1=1 ; IA2=IA1+NA-1 ; isA=true;
    Itotal=IA2;
end

if contains(CtrlVar.Inverse.InvertForField,'b')
    Ib1=Itotal+1 ; Ib2=Ib1+Nb-1 ; isb=true;
    Itotal=Ib2;
end

if contains(CtrlVar.Inverse.InvertForField,'B')
    IB1=Itotal+1 ; IB2=IB1+NB-1 ; isB=true;
    Itotal=IB2;
end

if contains(CtrlVar.Inverse.InvertForField,'C')
    IC1=Itotal+1 ; IC2=IC1+NC-1 ; isC=true;
    Itotal=IC2;
end

%
%   p = log(f)   <=> f=10^p
%
%   or
%
%   f=M^{1/2) p
%


if isA
    if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
        F.AGlen=10.^p(IA1:IA2);
    else
        F.AGlen=p(IA1:IA2);
    end
end

 if isb
     error('fdsa')
%     F.h=F.s-p(Ib1:Ib2) ;
%     F.B=p.*F.GF.node+(1-F.GF.node).*F.B;
%     
%     
%     %  bfloat=F.S - F.rho.*(F.s-p) /F.rhow;
%     %  dbfloat/dp= F.rho./F.rhow
%     %
%     % F.h=F.GF.node.*(F.s-F.B)+(1-F.GF.node).*F.hInit ;
%     
%     [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,[],F.h,F.S,F.B,F.rho,F.rhow);
 end

if isB
    
    % tested and works express geometrical variables in terms of p
    
    % here the control variable is B so I need to express all other variables in
    % terms of B this includes changing the thickness.
    %
    % The idea is to use s as given by measurements,i.e. s=Meas.s This gives b over
    % the floating areas from flotation argument using GF, and with it h and B Over
    % the grounded areas, again use s=Meas.s, and set b=B. However, this raises the
    % possibility that:
    %   1) As b is modified over the floating areas we get b<B in places. 2) As B is
    %   modified over the grounding areas, some areas go afloat.
    % 
    % This is a difficult one, but I assume Meas.s is know quite accuratly so
    % whatever is done, keep s=Meas.s. Also do not change B as it is given by p.
    %
    % So I must calculate b from s over the floating areas, and set b=B over the
    % grounded areas. However, this might cause b to be below B (b<B) for two
    % separate reasons: 1) Calculating b from the floating condition might results
    % in b<B. And changes in B below the GL (regularisation effects) may also give
    % rise to b<B.
    %
    % Solution:
    %
    %   1) Initially keeping GF, set b=B over grounded
    %      areas. But note that this will violate the GF mask if B is lowered as
    %      areas around the GL go afloat. If B is raised where grounded then locally
    %      GF does not change.
    %
    %   2) set b=max(b,B). Because of step 1, this will not change b any further
    %   over the grounded areas. This will also not change b over the floating areas
    %   provided B does not change there. But B may change over the grounded areas
    %   due to regularisation, and therefore this step is required. Therefore, b may
    %   be shifted upwards over the areas previously afloat and this will cause new
    %   areas to become grounded.
    %      
    %   3) set h=s-b
    %
    %   4) to ensure consistency with floating condition, recalculate GF from S and
    %   B for the new h.  Then recalculate b from Meas.s where afloat acording to
    %   the new GF, and then finally set again h=s-b.
    % 
    %   If B  is lowered over (previously) grounded areas, and therefore b as well, causing
    %   grounded areas to go afloat based on the new increased thickness h=Meas.s-b,
    %   then the recalculation of b in step 4 will only raise b, and hence not give
    %   rise to b<B situations. 
    %
    %   If B is lowered over previously floating areas, GF is not affected and there
    %   are no changes to b in step 4.
    %
    %   If B is shifted upwards over previously grounded areas, then GF 
    %   also does not change and there is no change in b in step 4. 
    %
    %   If B is shifted upwards over previous floating areas, causing possible
    %   grounding, then the recalculation of GF in step 4 ensures that b is only
    %   recalculated over areas previously afloat.
    %
    
    
    F.B=p(IB1:IB2) ;
    
    %F.B=F.GF.node.*p(IB1:IB2)+(1-F.GF.node).*Priors.B ;
    if CtrlVar.Inverse.OnlyModifyBedUpstreamOfGL
        [F.GF,GLgeo,GLnodes,GLele]=IceSheetIceShelves(CtrlVar,MUA,F.GF,GLgeo,GLnodes,GLele) ;
        F.B(~F.GF.NodesUpstreamOfGroundingLines)=Priors.B(~F.GF.NodesUpstreamOfGroundingLines) ;
    end
    
    F.s=Meas.s ; % note that since I'm not inverting for s, I must keep s fixed,
    % therefore calculate F.b over the floating areas from F.s using the floating relationship.
    
    [F.b,F.h,F.GF]=Calc_bh_From_sBS(CtrlVar,MUA,F.s,F.B,F.S,F.rho,F.rhow,F.GF); %
    
 
    
        
end

if isC

    if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
        F.C=10.^p(IC1:IC2);
    else
        F.C=p(IC1:IC2);
    end
    
end



end