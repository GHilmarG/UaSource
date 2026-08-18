
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
    fprintf("A sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tA)
    scale=log(10)*F.AGlen;
    dudA=dudA.*scale ; % using implicit expansion
    dvdA=dvdA.*scale ; % using implicit expansion
    if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)
        dhdA=dhdA.*scale ;
    else
        dhdA=[];
    end
end

if contains(CtrlVar.Inverse.InvertFor,"-B-")
    tB=tic;
    [dudB,dvdB,dhdB]=duvdBFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv) ;  % this has been tested against finite-differences and is good
    tB=toc(tB);
    fprintf("B sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tB)
end

if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
    tC=tic;
    [dudC,dvdC,dhdC]=duvdCFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv) ; % this has been tested against finite-differences and is good
    tC=toc(tC);
    fprintf("C sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tC)

    scale=log(10)*F.C;  % this has to be a row vector
    dudC=dudC.*scale ;   % using implicit expansion
    dvdC=dvdC.*scale ;   % using implicit expansion
    if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)
        dhdC=dhdC.*scale ;
    else
        dhdC=[];
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

    %% Test


    Funperturbed=F;

    NodeTest=randi(MUA.Nnodes);

    %% A
    if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)


   
        F=Funperturbed;
        DeltaRel=0.01;
        
        Field="AGlen";
        [dudApert,dvdApert,dhdotdApert]=FiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel);

        PlotModelAndFiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel,dudA,dvdA,dhdA,dudApert,dvdApert,dhdotdApert);


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
        PlotModelAndFiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel,dudB,dvdB,dhdB,dudBpert,dvdBpert,dhdotdBpert);


    end

    %% C
    if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)

    
        F=Funperturbed;
        DeltaRel=0.1;
        Field="C";
        [dudCpert,dvdCpert,dhdotdCpert]=FiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel);

        PlotModelAndFiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel,dudC,dvdC,dhdC,dudCpert,dvdCpert,dhdotdCpert);





        %%
    end

    drawnow
    input("duv_hdABC: Inspect, and then press RET to continue")

end