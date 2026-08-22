







function [dudField,dvdField,dhdotdField]=FiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,Node,DeltaRel,DeltaAbs)


if nargin< 11
    DeltaAbs=nan;
end

Field0=F.(Field);


% I do the perturbation in lin space, not in log space
%
%
if isnan(DeltaAbs)
    DeltaField=F.(Field)(Node)*DeltaRel;  % perturbation
else
    DeltaField=DeltaAbs;
end


% positive perturbation 


F.(Field)(Node)=F.(Field)(Node)+DeltaField ;
[UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
[UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
up=F.ub; vp=F.vb; dhdtp=dhdt;

% back to original
F.(Field)=Field0;


%Negative perturbation 

F.(Field)(Node)=F.(Field)(Node)-DeltaField ;



[UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
[UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
um=F.ub; vm=F.vb; dhdtm=dhdt;


% differences 
dudField=(up-um)/(2*DeltaField) ;
dvdField=(vp-vm)/(2*DeltaField) ;
dhdotdField=(dhdtp-dhdtm)/(2*DeltaField) ;

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

%

if Field~="B"
    % Here I do the log scaling
    scale=log(10)*Field0(Node);
    dudField=dudField.*scale ;
    dvdField=dvdField.*scale ;
    dhdotdField=dhdotdField.*scale ;



end


end