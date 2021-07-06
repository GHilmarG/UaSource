function CtrlVar=NrOfIntegrationPoints(CtrlVar)
    %%
    %  [nip,niph,BoundaryEdge]=NrOfIntegrationPoints(CtrlVar)
    %  defines nr of integration points for diagnostic (nip) and prognostic (nihp) equations
    %
    % Possible Nr of integration points: 1,3,4,6,7,9,12,13,16,19,28,37
    %
    % very important for the h and the implicit uvh cases to have one higher order of integration!!!!
    % nip must be above lin/quadradic for correct results! Why?
    % 
    % 3 node:
    %        nip=1 and niph=1 in semi-implicit is unstable
    %        nip=1 and niph=3 in semi-implicit is stable
    %        nip=3 and niph=3 in semi-implicit is stable
    % 6 node: 
    %        nip=4 in implicit uvh transient calculation leads to singular system
    %        nip=4 and niph=4  in semi-implicit transient calculation leads to singular system
    %        nip=4 and niph=6  in semi-implicit is fine
    %
    % 10 node;
    %        nip=7, niph=7 unstable in semi-implicit and implicit scemes 
    %        nip=7, niph=12 stable in semi-impicit scemes
    %        nip=12, niph=12 stable in both semi and implicit
    %
    %  The rule is that niph must be >=  3, 6, and 12 for 3, 6 and 10 node elements, respectively
    %  i.e. the prognostic equation must always be integrated with 1-order higher integration scheme
    % 
    %%

    persistent nCalls  nipLast niphLast
   


    if isempty(nCalls) ; nCalls=0 ; end
    
    nCalls=nCalls+1; 
    
    switch CtrlVar.TriNodes
        case 3 % minimum of 1 needed for a linear problem
            

            %nip=1 ; niph=1 ;  % not OK for both semi and fully implicit
            %nip=1; niph=3;    % OK for semi-implicit
            % nip=3; niph=3;   % OK for implicit and semi-implicit 
            
            %niph=3;
            %niph=4; % for some odd reason niph=4 sometimes causes convergence problems, while using niph=3 or niph=6 does not
            %nip=1 ; niph=3 ; % test
            %nip=3 ; niph=4 ; % test
            
            nip=6; niph=6;    % changed from 3 to 6 as default, 19 August, 2018
            
            
        case 6   % minimum of 4 needed for a linear problem
            
            %nip=4 ; % results in a singular system for implicit uvh
            %nip=6;
            %nip=6;   %  for GL problems overintegrating leads to impoved rates of convergence
            %nip=12;
            
            %nip=6 ; niph=6;
            % nip=7 ; niph=7;  
            
            %nip=9 ; niph=9;   % changed from 7 to 9 as default 19, August, 2018
            nip=12 ; niph=12;  % and then changed from 9 to 12 on 18 March, 2019. The reason for this increar, ie form 9 to 12 was that a few situations where found were 
                               % nip=9 was not sufficient for full
                               % convergence in a uv inversion. 
            
            
            
            
        case 10 % minimum of 7 needed for a linear problem
            % nip=7;
            % nip=7; niph=12;
            % nip=12 ; niph=12; 
             nip=16; niph=16;
            
            
        otherwise
            error(' case not recognised, TriNodes value incorrect')
    end
    
    %if CtrlVar.Implicituvh ; nip=niph ; end
    
   if ~isfield(CtrlVar,'nip') ; CtrlVar.nip=nip ; end
   if ~isfield(CtrlVar,'niph') ; CtrlVar.niph=niph ; end
    
    
    
    if isempty(CtrlVar.nip) ;  CtrlVar.nip=nip ; end
    if isempty(CtrlVar.niph) ;  CtrlVar.niph=niph ; end
   
    if nCalls==1 || nipLast~=nip || niphLast~=niph
        fprintf('Number of integration points: nip=%-i niph=%-i \n',CtrlVar.nip,CtrlVar.niph)
    end
    
    
    nipLast=nip ; niphLast=niph;
    
    
    
    
end