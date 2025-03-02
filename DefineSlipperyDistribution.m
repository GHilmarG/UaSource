
function  [UserVar,C,m,q,muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,F)
    
    
    %%
    %
    % [UserVar,C,m,q,muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,F)
    %
    % [UserVar,C,m,q,muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)
    %
    %
    % Defines sliding-law parameters.
    %
    % The sliding law used is determined by the value of 
    %
    %   CtrlVar.SlidingLaw
    %
    % which is defined in 
    %
    %   DefineInitialInputs.m
    %
    % See description in Ua2D_DefaultParameters.m for further details and the
    % UaCompendium.pdf.
    %
    %%
    
    
    m=3;
    C0=3.16e6^(-m)*1000^m*365.2422*24*60*60;
    
    C=C0+zeros(MUA.Nnodes,1);
    
    
    q=1 ;      % only needed for Budd sliding law
    muk=0.5 ;  % required for Coulomb friction type sliding law as well as Budd, minCW (Tsai), rCW  (Umbi) and rpCW (Cornford. 
    
    
end
