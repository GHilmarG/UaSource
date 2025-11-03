function lambdaStar=CalcLambdaStar(CtrlVar,MUA,L,lambda)

% Calculates 'physical' Lagrange parameter values    
    
% MLC=BCs2MLC(CtrlVar,MUA,BCs1) ; L=MLC.hL ; M=MUA.M ;   

M=MUA.M; 

lambdaStar=(L*L')\L*(M\(L'*lambda)) ;



end