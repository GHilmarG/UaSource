





function [ab,dabdh]=DraftDependentMeltParameterisations(UserVar,CtrlVar,F,MRP)

%%
%
%   [ab,dabdh]=DraftDependentMeltParameterisations(UserVar,CtrlVar,F,MRP)
%
%
% Returns basal melt rates, based on several depth-dependent parameterisations.
%
%
%  ab is :    abMin for F.b>dMin
%             abMax for F.b<dMax
%
%  and varies linearly in between.
%
%
% These types of parameterisations have been used in various papers. Labels might not correspond exactly to the labels used in
% previous publications.
%
%
% Note !!!
% The user will in addition need to make sure that melt is only applied over floating nodes, downstream of the grounding line
% This can be done in DefineMassBalance.m as follows:
%
%   [LakeNodes,OceanNodes]=LakeOrOcean3(CtrlVar,MUA,GF)
%   ab(~OceanNodes) = 0;
%   dabdh(~OceanNodes)=0;
%
%%


switch MRP


    case {"0","l0"}

        fprintf(' MeltRate0 \n ');
        dMin=-400 ;  abMin=0 ;
        dMax=-800 ;  abMax=-100  ;

    case {"1","l1"}

        fprintf(' MeltRate1 \n ');
        dMin=-400 ;  abMin=0 ;
        dMax=-800 ;  abMax=-200  ;


    case {"2","l2"}

        fprintf(' MeltRate2 \n ');
        dMin=-200 ;  abMin=0 ;
        dMax=-800 ;  abMax=-100  ;


    case {"3","l3"}

        fprintf(' MeltRate3 \n ');
        dMin=-200 ;  abMin=0 ;
        dMax=-800 ;  abMax=-200  ;

    case {"4","l4"}

        % fprintf(' MeltRate4 \n ');
        dMin=-400 ;  abMin=0 ;
        dMax=-500 ;  abMax=-50  ;

    case {"5","l5"}

        fprintf(' MeltRate5 \n ');
        dMin=-400 ;  abMin=0 ;
        dMax=-600 ;  abMax=-100  ;

    case {"ASE1","lASE1"}  % an atempt to come up with something for the Amunsen Sea Embayment that works for the whole region

        % Here the thermocline transition range is kept constant across the boundary
        % But the melt in the thermocline is dependend on the y ups values

        dMin=-400 ; abMin=0 ;
        dMax=-500  ; 

        % abMax(y=-400e3)=-100;
        % abMax(y=-700)=-50  ;
        
        % abMax=a0+a1* (y+700) ; % -> a0=-50
        % abMax=-50 + a1* (-400+700) = -50+a1 *300 = -100 ; -> a1=-50/300e3 ;

        abMax=-50-(50/300e3)*(700e3+F.y) ; 



    otherwise

        fprintf("DraftDependentMeltParameterisations : case not found")
        dMin=nan ;  abMin=nan ;
        dMax=nan ;  abMax=nan  ;


end



F.as=zeros(size(F.h)) ;
F.ab=zeros(size(F.h)) ;
F.dasdh=zeros(size(F.h)) ;
F.dabdh=zeros(size(F.h)) ;


if isscalar(abMax)
    abMax=abMax+F.x*0;
end



if contains(MRP,"l")  % this is a smoothed version of the melt paramerisations

    b0=(dMin+dMax)/2 ;
    k=2/(dMin-dMax) ;
    ab=abMin+    (1-HeavisideApprox(k,F.b,b0)).*(abMax-abMin) ;

    dabdh=-DiracDelta(k,F.b,b0).*(abMax-abMin).*(-F.rho/F.rhow) ;
    %   dab/db  db/dh

else

    dabdh=zeros(size(F.x)) ;
    ab=zeros(size(F.x)) ;


    % Note, the dab/dh calculation is not exact, as it is missing a Dirac delta
    % term

    I=F.b>dMin ; ab(I)=abMin; dabdh(I)=0;        % above dMin

    I= F.b<= dMin & F.b >= dMax ;

    % b= -h rho/rhow
    ab(I)=abMax(I).*(F.b(I)-dMin)/(dMax-dMin); % negative values, because it is a basal ablation
    dabdh(I)=(-F.rho(I)/F.rhow)  .* (abMax(I)/(dMax-dMin));

    I=F.b<dMax ; ab(I)=abMax(I); dabdh(I)=0;        % below dMax


end


% Note !!!
% The user will in addition need to make sure that melt is only applied over floating nodes, downstream of the grounding line
% This can be done in DefineMassBalance.m as follows:
%
%   [LakeNodes,OceanNodes]=LakeOrOcean3(CtrlVar,MUA,GF)
%   ab(~OceanNodes)=0;
%   dabdh(~OceanNodes)=0;
%
%%


end