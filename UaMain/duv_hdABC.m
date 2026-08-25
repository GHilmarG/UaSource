
function [dudA,dvdA,dhdA,dudB,dvdB,dhdB,dudC,dvdC,dhdC]=duv_hdABC(UserVar,CtrlVar,RunInfo,MUA,F,l,BCs)


narginchk(7,7)
nargoutchk(9,9)

dudA=[]; dvdA=[]; dhdA=[];
dudB=[]; dvdB=[]; dhdB=[];
dudC=[]; dvdC=[]; dhdC=[];

[UserVar,RunInfo,F,l,KdFuvduv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);

% CtrlVar.uvAssembly.ZeroFields=false;
% CtrlVar.uvMatrixAssembly.Ronly=false;
% [~,KdFuvduv]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);


if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)
    tA=tic;
    [dudA,dvdA,dhdA]=duvhdotdAFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv) ;  % this has been tested against finite-differences and is good, also for dhdotdA
    tA=toc(tA);
    %fprintf("A sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tA)
 
    ln10 = log(10);
    ScaleMatrix=spdiags(F.AGlen(:)*ln10, 0, MUA.Nnodes, MUA.Nnodes);
    dudA = dudA * ScaleMatrix ; 
    dvdA = dvdA * ScaleMatrix ; 

    if ~isempty(dhdA)
        dhdA=dhdA*ScaleMatrix ;
    end
end

if contains(CtrlVar.Inverse.InvertFor,"-B-")
    tB=tic;
    [dudB,dvdB,dhdB]=duvdBFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv) ;  % this has been tested against finite-differences and is good
    tB=toc(tB);
    %fprintf("B sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tB)
end

if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
    tC=tic;
    [dudC,dvdC,dhdC]=duvdCFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv) ; % this has been tested against finite-differences and is good
    tC=toc(tC);
    %fprintf("C sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tC)


    ln10 = log(10);
    ScaleMatrix=spdiags(F.C(:)*ln10, 0, MUA.Nnodes, MUA.Nnodes);
    dudC = dudC * ScaleMatrix ; 
    dvdC = dvdC * ScaleMatrix ; 

    if ~isempty(dhdC)
        dhdC=dhdC*ScaleMatrix ;
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

    FiniteDifferenceTestAndPlots(F,MUA,CtrlVar,UserVar,RunInfo,BCs,l,dudA,dvdA,dhdA,dudB,dvdB,dhdB,dudC,dvdC,dhdC);

end

function FiniteDifferenceTestAndPlots(F,MUA,CtrlVar,UserVar,RunInfo,BCs,l,dudA,dvdA,dhdA,dudB,dvdB,dhdB,dudC,dvdC,dhdC)
%% Test


Funperturbed=F;

NodeTest=randi(MUA.Nnodes);

%% A
if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)



    F=Funperturbed;
    DeltaRel=1e-4;

    Field="AGlen";
    [dudApert,dvdApert,dhdotdApert]=FiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel);

    SubtitleString="sensitivites are with respect to $\log_{10}A$";

    PlotModelAndFiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel,dudA,dvdA,dhdA,dudApert,dvdApert,dhdotdApert,SubtitleString);

end


%% B
if contains(CtrlVar.Inverse.InvertFor,"-B-")
    %% du/dB
    F=Funperturbed;

    CtrlVar.Calculate.Geometry="bh-FROM-sBS" ;


    F=Funperturbed;
    DeltaRel=nan;
    DeltaAbs=1;
    Field="B";
    [dudBpert,dvdBpert,dhdotdBpert]=FiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel,DeltaAbs);

    % figures
    Field="B";
    SubtitleString="sensitivites are with respect to $B$";

    PlotModelAndFiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel,dudB,dvdB,dhdB,dudBpert,dvdBpert,dhdotdBpert,SubtitleString);

end

%% C
if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)


    F=Funperturbed;
    DeltaRel=1e-4;
    Field="C";
    [dudCpert,dvdCpert,dhdotdCpert]=FiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel);

    SubtitleString="sensitivites are with respect to $\log_{10}C$"; 

    PlotModelAndFiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel,dudC,dvdC,dhdC,dudCpert,dvdCpert,dhdotdCpert,SubtitleString);

    %%
end

