





function [KdudA,KdvdA,KdudB,KdvdB,KdudC,KdvdC]=duv_hdABC(CtrlVar,MUA,F,l,BCs,isA,isB,isC)

%% This returns Jacobian, i.e. du/dA and so on.
%
%
%
% The mesh-independent object is the kernel $K(x,y)$
%
% $$\delta u(x)= \int K(x,y) \, \delta B(y) \; dA_y $$
%
%
% $$\frac{\partial u_i}{\partial B_j} = \int K(x,y) \, \phi_j(y) \; dA_y $$
%
% or
%
% $$K = J M^{-1} $$
%
% where $K$ is the matrix kernel and $J$ the Jacobian and $M$ the mass matrix excellent
%
% Units check:  For $B$: $J$ is $yr^{-1}$ $K$ is $yr^{-1}\,m^{-2}$ and $\int K \, \delta B \, dA $ give $m/yr$.
%
%
%%


narginchk(5,8)
nargoutchk(6,6)


KdudA=[]; KdvdA=[]; KdhdA=[];
KdudB=[]; KdvdB=[]; KdhdB=[];
KdudC=[]; KdvdC=[]; KdhdC=[];

[~,~,F,l,KdFuvduv]= uv([],[],CtrlVar,MUA,BCs,F,l);

if nargin<=5

    [isA,isB,isC] = isABC(CtrlVar);

end

if isA
    %tA=tic;
    [KdudA,KdvdA,KdhdA]=duvhdotdAFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv) ;  % this has been tested against finite-differences and is good, also for dhdotdA
    %tA=toc(tA);
    %fprintf("A sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tA)

    ln10 = log(10);
    ScaleMatrix=spdiags(F.AGlen(:)*ln10, 0, MUA.Nnodes, MUA.Nnodes);
    KdudA = KdudA * ScaleMatrix ;
    KdvdA = KdvdA * ScaleMatrix ;

    if ~isempty(KdhdA)
        KdhdA=KdhdA*ScaleMatrix ;
    end
end

if isB
    %tB=tic;
    % [dudB,dvdB,dhdB]=duvdBFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv) ;  % this has been tested against finite-differences and is good
    [KdudB,KdvdB]=duvdBGeneral(CtrlVar,MUA,F,l,BCs,KdFuvduv) ; % this is the general case, not assuming the ice to be grounded everywhere
    %tB=toc(tB);
    % fprintf("B sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tB)
end

if isC
    %tC=tic;
    [KdudC,KdvdC,KdhdC]=duvdCFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv) ; % this has been tested against finite-differences and is good
    %tC=toc(tC);
    %fprintf("C sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tC)


    ln10 = log(10);
    ScaleMatrix=spdiags(F.C(:)*ln10, 0, MUA.Nnodes, MUA.Nnodes);
    KdudC = KdudC * ScaleMatrix ;
    KdvdC = KdvdC * ScaleMatrix ;

    if ~isempty(KdhdC)
        KdhdC=KdhdC*ScaleMatrix ;
    end
end




%%
% log10 sensitivities
%
% du/dA=du/dx  dx/dA
%
% x=ln(A)
%
% du/dA=du/dx  d(ln(A))/dA  = du/dx   1/A
%
% Therefore
%
% du/d(ln(A)) = A du/dA
%
% or
%
% du/d(ln(A)) = log(10) A du/dA
%
%


if CtrlVar.Inverse.TestDirectAdjoint.isTrue

    % Note: This works fine for all perturbations. However, if the idea is to test (du/dB,dv/dB) a better bespoke test is
    % implemented in TestduvdBGeneral.m
    FiniteDifferenceTestAndPlots(F,MUA,CtrlVar,BCs,l,KdudA,KdvdA,KdhdA,KdudB,KdvdB,KdhdB,KdudC,KdvdC,KdhdC);

end

function FiniteDifferenceTestAndPlots(F,MUA,CtrlVar,BCs,l,dudA,dvdA,dhdA,dudB,dvdB,dhdB,dudC,dvdC,dhdC)
%% Test


Funperturbed=F;

NodeTest=randi(MUA.Nnodes);

%% A
if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)



    F=Funperturbed;
    DeltaRel=1e-4;

    Field="AGlen";
    [dudApert,dvdApert,dhdotdApert]=FiniteDifferenceSensitivities(CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel);

    SubtitleString="sensitivites are with respect to $\log_{10}A$";

    PlotModelAndFiniteDifferenceSensitivities(CtrlVar,MUA,BCs,F,l,Field,NodeTest,dudA,dvdA,dhdA,dudApert,dvdApert,dhdotdApert,SubtitleString);

end


%% B
if contains(CtrlVar.Inverse.InvertFor,"-B-")
    %% du/dB



    F=Funperturbed;
    DeltaRel=nan;
    DeltaAbs=0.1;
    Field="B";
    [dudBpert,dvdBpert,dhdotdBpert]=FiniteDifferenceSensitivities(CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel,DeltaAbs);

    % figures
    Field="B";
    SubtitleString="sensitivites are with respect to $B$";

    PlotModelAndFiniteDifferenceSensitivities(CtrlVar,MUA,BCs,F,l,Field,NodeTest,dudB,dvdB,dhdB,dudBpert,dvdBpert,dhdotdBpert,SubtitleString);
%%
end

%% C
if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)


    F=Funperturbed;
    DeltaRel=1e-2;
    Field="C";
    [dudCpert,dvdCpert,dhdotdCpert]=FiniteDifferenceSensitivities(CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel);

    SubtitleString="sensitivites are with respect to $\log_{10}C$";

    PlotModelAndFiniteDifferenceSensitivities(CtrlVar,MUA,BCs,F,l,Field,NodeTest,dudC,dvdC,dhdC,dudCpert,dvdCpert,dhdotdCpert,SubtitleString);

    %%
end

