function [UserVar,RunInfo,F,F0,l,Kuv,Ruv,Lubvb]= WTSHTF(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,Fm1,l)

%%
%
% If the implicit uvh solution did not converge, take a semi-implicit step and
% continue.
%
%  n adapt time step the dt will be reduced by the fraction
%
%     CtrlVar.ATStimeStepFactorDownNOuvhConvergence
%
% So this should hopefully cause a graceful reduction in time step without
% significant reduction in accuracy.
%


filename="Dumpfile_Ua2D-"+CtrlVar.Experiment+".mat";
fprintf(' ===>>> uvh did not converge! Saving all data in a dumpfile %s \n',filename)
save(filename)

% Make sure not to overwrite the RunInfo from the uvh step

% Since 


% [UserVar,~,F,F0,l,Kuv,Ruv,Lubvb]= uvhSemiImplicit(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,Fm1,l);

[UserVar,RunInfo,F,F0,l,Kuv,Ruv,Lubvb,duv1NormVector]= uvhSemiImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F0,l,BCs)  ; 

end


