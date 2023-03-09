function [tau,ECN,K]=SUPGtau(CtrlVar,v,l,dt,tauOption,multiplier)


%
% v    : speed
%


if nargin < 6
    multiplier=1;
end

taut=dt/2+eps;   % temporal definition
taus=l./(2*v+CtrlVar.SpeedZero);  % spatial definition
% ECN=v.*dt./l;

ECN=taut./taus ;


%%
%
% tau1 : often recomended in textbooks for linear diffusion equations with
%        spatially constant non-zero advection velocity
% taut : dt/2,  'temporal' definition, independed of velocity
% taus : l/(2u) 'spatial definition', independent of time step
% tau2 : 1./(1./taut+1./taus), an 'inverse' average of taus and taut


% And now I must consider the possibility that speed is zero, in which case
% the above expression fails and must be replaced by the correct limit which is
% tau1 -> dt/6 as speed -> 0


switch tauOption

    case "tau1"   %  typical textbook recomendation for spatially constant (and non-zero) speed for linear advection equation

        K=coth(ECN)-1./ECN;  % (1/ECN+ECN/3+..) -1/ECN=ECN/3  if ECN->0
        % turns out the expression for K starts to suffer from numerical errors for ECN < 1e-6
        % x=logspace(-10,-5); figure ; semilogx(x,coth(x)-1./x,'r') ; hold on ; semilogx(x, x/3,'b')
        I=ECN < 1e-6 ; K(I)=ECN(I)/3 ;  % replaced by the Taylor expansion
        tau=K.* taus ; % l./v/2;
        %I=v<100*eps ; tau(I)=dt/6;

    case "tau2"   %  inversly weighted average of spatial and temporal tau
        % tau=1./(1./taut+1./taus);
        %    =1./(1./(dt/2)  + 1/(l/2v))
        tau=(dt/2).*1./(1+ECN) ;

    case "taus"   % 'spatial' definition, independent of time step
        tau=taus ;

    case "taut"   % 'temporal' definition, independent of speed
        tau=taut;

    case "tau0"   % set to zero...!

        tau=zeros(size(v),'like',v);

    otherwise
        error("in SUPGtau case not found")
end



tau=multiplier*tau ;

if CtrlVar.doplots  && CtrlVar.PlotSUPGparameter && ~isempty(MUA)

    
    K=coth(ECN)-1./ECN;  % (1/ECN+ECN/3+..) -1/ECN=ECN/3  if ECN->0
    I=ECN < 1e-6 ; K(I)=ECN(I)/3 ;  % replaced by the Taylor expansion
    tau1=K.* taus ; % l./v/2;

    tau2=(dt/2).*1./(1+ECN) ;
    taut=taut+zeros(size(v),'like',v);

    FindOrCreateFigure("taus SUPG") ; histogram(taus) ; hold on ; histogram(taut) ; histogram(tau1) ; histogram(tau2) ; legend

    % incomplete because MUA not handed over in the integration point loop
%     FindOrCreateFigure("taus SUPG") ;
%     subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,taut) ; title('taut')
%     subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,taus) ; title('taus')
%     subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,tau1) ; title('tau1')
%     subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,tau2) ; title('tau2')

end





end