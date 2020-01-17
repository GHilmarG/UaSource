function [UserVar,RunInfo,r,ruv,rh,rl,R,K,frhs,grhs]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,Fext0)



narginchk(15,15)

Nout=nargout;

F1.ub=F1.ub+gamma*dub;
F1.vb=F1.vb+gamma*dvb;
F1.h=F1.h+gamma*dh;
luvh=luvh+gamma*dl;

CtrlVar.uvhMatrixAssembly.ZeroFields=false;
if Nout<8
    CtrlVar.uvhMatrixAssembly.Ronly=true;
else
    CtrlVar.uvhMatrixAssembly.Ronly=false;
end

[UserVar,RunInfo,R,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1);


if ~isempty(L)
    
    frhs=-R-L'*luvh;
    grhs=cuvh-L*[F1.ub;F1.vb;F1.h];
    
    %frhs=-R-L'*(luv+gamma*dl);
    %grhs=cuvh-L*[ub+gamma*dub;vb+gamma*dvb;h+gamma*dh];
else
    frhs=-R;
    grhs=[];
end

% Testing:
%
% in a uv solution the uv residuals have different physical dimentions from h.
% One way of dealing with this is to multiply in a transient simulation uv with
% dt, or divide h with dt
%
% Nuv=2*MUA.Nnodes; frhs(Nuv+1:end)=frhs(Nuv+1:end)/CtrlVar.dt;

[r,rl,ruv,rh]=ResidualCostFunction(CtrlVar,MUA,L,frhs,grhs,Fext0,"-uvh-");



if isequal(CtrlVar.Residual.uvh,'uv')
    r=ruv;
end

end

