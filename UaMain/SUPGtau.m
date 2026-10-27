





function [tau,ECN,K]=SUPGtau(CtrlVar,v,l,dt,tauOption,multiplier)


%
% v    : speed
% l    : element size  
%
% ECN : element Peclet number, ratio of distance traveled, dt v, to the
% element size, l. 
%
% tauOption : "tau1", "tau2", "taus", "taut" or "tau0" (CtrlVar.uvh.SUPG.tau), or the corresponding numeric code 1-5, see
%             SUPGtauCode.m (eg determined once per assembly in the compiled assembly, UaSource/UaCompiled)
%
% K : only defined for "tau1" (NaN otherwise)
%
% (26 Oct 2026) Written so that MATLAB Coder can use it directly (the compiled uvh assembly): the option can be given as a numeric
% code, and the plotting is excluded from compiled code.


if nargin < 6
    multiplier=1;
end

% the dimension of tau is time

taut=dt/2+eps;   % temporal definition
taus=l./(2*v+CtrlVar.SpeedZero);  % spatial definition
% ECN=v.*dt./l;

ECN=taut./taus ; % ECN=taut/taus=(dt/2)/(l/(2 v)) = dt v / l  


%%
%
% tau1 : often recommended in textbooks for linear diffusion equations with
%        spatially constant non-zero advection velocity
% taut : dt/2,  'temporal' definition, independent of velocity
% taus : l/(2u) 'spatial definition', independent of time step
% tau2 : 1./(1./taut+1./taus), an 'inverse' average of taus and taut


% And now I must consider the possibility that speed is zero, in which case
% the above expression fails and must be replaced by the correct limit which is
% tau1 -> dt/6 as speed -> 0


if isnumeric(tauOption)
    Code=tauOption ;
else
    Code=SUPGtauCode(tauOption) ;
end

switch Code

    case 1   % "tau1"   %  typical textbook recommendation for spatially constant (and non-zero) speed for linear advection equation

        K=coth(ECN)-1./ECN;  % (1/ECN+ECN/3+..) -1/ECN=ECN/3  if ECN->0
        % turns out the expression for K starts to suffer from numerical errors for ECN < 1e-6
        % x=logspace(-10,-5); figure ; semilogx(x,coth(x)-1./x,'r') ; hold on ; semilogx(x, x/3,'b')
        I=ECN < 1e-6 ; K(I)=ECN(I)/3 ;  % replaced by the Taylor expansion
        tau=K.* taus ; % l./v/2;
        %I=v<100*eps ; tau(I)=dt/6;

    case 2   % "tau2"   %  inversely weighted average of spatial and temporal tau
        K=NaN(size(ECN),'like',ECN);
        % tau=1./(1./taut+1./taus);
        %    =1./(1./(dt/2)  + 1/(l/2v))
        tau=(dt/2).*1./(1+ECN) ;

    case 3   % "taus"   % 'spatial' definition, independent of time step
        K=NaN(size(ECN),'like',ECN);
        tau=taus ;

    case 4   % "taut"   % 'temporal' definition, independent of speed
        K=NaN(size(ECN),'like',ECN);
        tau=taut;

    case 5   % "tau0"   % set to zero...!
        K=NaN(size(ECN),'like',ECN);

        tau=zeros(size(v),'like',v);

    otherwise
        error("in SUPGtau case not found")
end



tau=multiplier*tau ;

% (26 Oct 2026) not in compiled code; the condition ~isempty(MUA) removed (MUA is not an input, so this failed, and it is not needed
% for the histograms)
if coder.target('MATLAB') && CtrlVar.doplots  && CtrlVar.PlotSUPGparameter

    
    K=coth(ECN)-1./ECN;  % (1/ECN+ECN/3+..) -1/ECN=ECN/3  if ECN->0
    I=ECN < 1e-6 ; K(I)=ECN(I)/3 ;  % replaced by the Taylor expansion
    tau1=K.* taus ; % l./v/2;

    tau2=(dt/2).*1./(1+ECN) ;
    taut=taut+zeros(size(v),'like',v);

    fHist=FindOrCreateFigure("taus SUPG") ; clf(fHist)
    histogram(taus) ; 
    hold on ; 
    histogram(taut)
    histogram(tau1) ; 
    histogram(tau2) ; 
    set(gca, "XScale", "log")
    legend

    % incomplete because MUA not handed over in the integration point loop
%     FindOrCreateFigure("taus SUPG") ;
%     subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,taut) ; title('taut')
%     subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,taus) ; title('taus')
%     subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,tau1) ; title('tau1')
%     subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,tau2) ; title('tau2')

end





end