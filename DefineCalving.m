function [UserVar,LSF,c]=DefineCalving(UserVar,CtrlVar,MUA,F,BCs)
    
    
    %%
    %
    %   [UserVar,LSF,CalvingRate]=DefineCalving(UserVar,CtrlVar,MUA,F,BCs)
    %
    % Define calving the Level-Set Field (LSF) and the Calving Rate Field (c)
    %
    % Both the Level-Set Field (LSF) and the Calving-Rate Field (c) must be defined over
    % the whole computational domain.
    % 
    %
    % The LSF should, in general, only be defined in the beginning of the run and set the
    % initial value for the LSF. However, if required, the user can change LSF at any time
    % step. The LSF is evolved by solving the Level-Set equation, so any changes done to
    % LSF in this m-file will overwrite/replace the previously calculated values for LSF.
    %
    % The calving-rate field, c, is an input field to the Level-Set equation and needs to
    % be defined in this m-file in each call.
    %
    % The variable F has F.LSF and F.c as subfields. In a transient run, these will be the
    % corresponding values from the previous time step.
    % 
    % If you do not want to modify LSF,  set
    %
    %   LSF=F.LSF
    %
    %
    % Also, if you do not want to modify c, you could in prinicple set
    %
    %   c=F.c
    %
    % However, note that in contrast to LSF, c is never evolved by Úa.  (Think of c as an
    % input variable similar to the input as and ab for upper and lower surface balance,
    % etc.)
    %
    % Initilizing the LSF is the task of the user and needs to be done in this m-file.
    % Typically LSF is defined as a signed distance function from the initial calving
    % front position. There are various ways of doing this and you might find the matlab
    % function
    % 
    %   pdist2
    %
    % usefull to do this. Also look at
    %
    %   ReinitializeLevelSet.m
    %
    % for ideas on how to initialize the level set.
    %
    %%
    
    
    LSF=F.LSF ;  % always do this if you do not want to modify the LSF
    
    if CtrlVar.CurrentRunStepNumber< 2  % initialize the Level-Set-Field
        
        % LSF set equal to signed distance from x=xc
        xc=200e3;  % this is the initial calving front
        LSF=xc-MUA.coordinates(:,1) ;
        
    end
    
    q=-2;
    k=86322275.9814533 ;
    
    switch UserVar.Calving
        
        case "Function of analytical thickness"
            % First for testing, define calving rate as a function of the analytical thickness
            % profile.
            
            [s,b]=AnalyticalOneDimentionalIceShelf(CtrlVar,MUA);
            h=s-b;
            
            
            c=k*h.^q;
            
        case "Function of numerical thickness"
            
            c=k*F.h.^q;
            % The issue is that this goes to infinity with h to zero. 
            % 
            c(c>2000)=2000; 
    end
    
    c=c(:) ;
    c=c.*(1-F.GF.node) ;

end
