function [tau,tau1,tau2,taus,taut,ECN,K,l]=CalcSUPGtau(CtrlVar,EleAreas,u,v,dt,MUA)


%
%  Calculates nodal based tau values to be used in the SUPG method.
%
% Area=l^2/2 -> l=sqrt( 2 Area) 

narginchk(5,6)

if nargin==5
    MUA=[];
end


if numel(u)~=numel(EleAreas)  % this could be called over nodes or elements (i.e. integration points)
    %    error('afsd')  % ; I've got rid of almost these cases, but still used in dhdtExplicitSUPG
    % Here MUA must be given as an input
    if ~nargin==6
        error("CalcSUPGtau:IncorrectNumberOfInputArguments","Incorrect number of input arguements")
    end
    
    M=Ele2Nodes(MUA.connectivity,MUA.Nnodes);  % this is for a call over nodes, try to get rid of this use
    l=M*sqrt(2*MUA.EleAreas) ;
else
    l=sqrt(2*EleAreas) ;
end

speed=sqrt(u.*u+v.*v);

ECN=speed.*dt./l;  % non-dimentional

K=coth(ECN)-1./ECN;  % (1/ECN+ECN/3+..) -1/ECN=ECN/3  if ECN->0
% turns out the expression for K starts to suffer from numerical errors for ECN < 1e-6
% x=logspace(-10,-5); figure ; semilogx(x,coth(x)-1./x,'r') ; hold on ; semilogx(x, x/3,'b')
I=ECN < 1e-6 ; K(I)=ECN(I)/3 ;  % replaced by the Taylor expansion

%%
%
% tau1 : often recomended in textbooks for linear diffusion equations with
%        spatially constant non-zero advection velocity
% taut : dt/2,  'temporal' definition, independed of velocity
% taus : l/(2u) 'spatial definition', independent of time step
% tau2 : 1./(1./taut+1./taus), an 'inverse' average of taus and taut
tau1=K.*l./speed/2;

% And now I must consider the possibility that speed is zero, in which case
% the above expression fails and must be replaced by the correct limit which is
% tau1 -> dt/6 as speed -> 0
I=speed<100*eps ; tau1(I)=dt/6;

taut=dt/2+eps++zeros(size(u),'like',u);
taus=0.5*l./(speed+CtrlVar.SpeedZero);  % Now this must go down to zero gracefully...
tau2=1./(1./taut+1./taus);

switch CtrlVar.Tracer.SUPG.tau
    
    case 'tau1'   %  typical textbook recomendation for spatially constant (and non-zero) speed for linear advection equation
        tau=tau1;
    case 'tau2'   %  inversly weighted average of spatial and temporal tau
        tau=tau2;
    case 'taus'   % 'spatial' definition, independent of time step
        tau=taus;
    case 'taut'   % 'temporal' definition, indepenent of speed
        tau=taut ; 
    otherwise
        error('in CalcSUPGtau case not found')
end


if CtrlVar.doplots  && CtrlVar.PlotSUPGparameter && ~isempty(MUA)


     
    FindOrCreateFigure("taus SUPG") ;
    subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,taut) ; title('taut')
    subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,taus) ; title('taus')
    subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,tau1) ; title('tau1')
    subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,tau2) ; title('tau2')
end


%fprintf('tau1=%f \t tau2=%f \t taus=%f \t taut=%f \n ',mean(tau1),mean(tau2),mean(taus),mean(taut))


%%

end