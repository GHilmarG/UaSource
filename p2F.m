function F=p2F(CtrlVar,MUA,p,F)

narginchk(4,4)

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
if CtrlVar.Inverse.pPreMultiplier=="M"
    Area=TriAreaTotalFE(MUA.coordinates,MUA.connectivity);
    
    if isA ;  p(IA1:IA2)=Area*(MUA.M\p(IA1:IA2)); end
    if isb ;  p(Ib1:Ib2)=Area*(MUA.M\p(Ib1:Ib2)); end
    if isB ;  p(IB1:IB2)=Area*(MUA.M\p(IB1:IB2)); end
    if isC ;  p(IC1:IC2)=Area*(MUA.M\p(IC1:IC2)); end
    
end


if isA
    if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
        F.AGlen=10.^p(IA1:IA2);
    else
        F.AGlen=p(IA1:IA2);
    end
end

if isb
    F.h=F.s-p(Ib1:Ib2) ;
    F.B=p.*F.GF.node+(1-F.GF.node).*F.BInit;
    
    
    %  bfloat=F.S - F.rho.*(F.s-p) /F.rhow;
    %  dbfloat/dp= F.rho./F.rhow
    %
    % F.h=F.GF.node.*(F.s-F.B)+(1-F.GF.node).*F.hInit ;
    
    [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,[],F.h,F.S,F.B,F.rho,F.rhow);
end

if isB
    
    % tested and works
    % express geometrical variables in terms of p
    
    F.B=p(IB1:IB2) ;
    F.h= F.hInit.*(1-F.GF.node)  + F.GF.node.* (F.sInit - p)  ;   % h = s - b
    
    %         bfloat=F.S - F.rho.*F.h /F.rhow;
    %
    %         F.b=F.GF.node.*p + (1-F.GF.node) .* bfloat ;
    %            =F.GF.node.*p + (1-F.GF.node) .* (F.S - F.rho.*F.h /F.rhow)
    %            =F.GF.node.*p + (1-F.GF.node) .* (F.S - F.rho.*(F.hInit.*(1-F.GF.node)  + F.GF.node.* (F.sInit - p))  ;
    F.b=F.GF.node.*p + (1-F.GF.node) .* (F.S - F.rho.*( F.hInit.*(1-F.GF.node)  + F.GF.node.* (F.sInit - p)   )./F.rhow);
    
    % b = GF  p + (1-GF) (S-rho (h0 (1-GF) + GF (s0-p) /rhow)
    % b = GF  p + (1-GF) (S-rho (h0 (1-GF) + GF (s0-p) /rhow)
    % db/dp = GF + (1-GF) rho GF/rhow
    
    %         dB/dp = 1
    %         dh/dp = -GF.node
    %         db/dp = GF.node + (1-GF.node) (0 - rho dhdp/rhow)
    %               = GF.node + (1-GF.node) (0 -- rho GF.node/rhow)
    %               = GF.node + (1-GF.node) (+ rho GF.node/rhow)
    %               = GF.node + rho (1-GF.node) GF.node/ rhow
    %
    
    [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,[],F.h,F.S,F.B,F.rho,F.rhow);
    
    % F.h= F.hInit.*(1-F.GF.node)  + F.GF.node.* (F.sInit - p)  ;  % because GF has changed
    
    
end

if isC

    if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
        F.C=10.^p(IC1:IC2);
    else
        F.C=p(IC1:IC2);
    end
    
end



end