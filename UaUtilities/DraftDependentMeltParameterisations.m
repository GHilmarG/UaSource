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
%  and varies linearly inbetween.
%
%
% These types of parameterisations have been used in various papers. Labels might not correspond exactly to the labels used in
% previous publicatioins.
% 
%
%%

switch MRP


    case "0"

        fprintf(' MeltRate0 \n ');
        dMin=-400 ;  abMin=0 ;
        dMax=-800 ;  abMax=-100  ;

    case "1"

        fprintf(' MeltRate1 \n ');
        dMin=-400 ;  abMin=0 ;
        dMax=-800 ;  abMax=-200  ;


    case "2"

        fprintf(' MeltRate2 \n ');
        dMin=-200 ;  abMin=0 ;
        dMax=-800 ;  abMax=-100  ;


    case "3"

        fprintf(' MeltRate3 \n ');
        dMin=-200 ;  abMin=0 ;
        dMax=-800 ;  abMax=-200  ;

    case "4"

        fprintf(' MeltRate4 \n ');
        dMin=-400 ;  abMin=0 ;
        dMax=-500 ;  abMax=-50  ;

    case "5"

        fprintf(' MeltRate5 \n ');
        dMin=-400 ;  abMin=0 ;
        dMax=-600 ;  abMax=-100  ;

    otherwise

        fprintf("DraftDependentMeltParameterisations : case not found")
        dMin=nan ;  abMin=nan ;
        dMax=nan ;  abMax=nan  ;


end


if ~isequal(F.as,F.x)
    F.as=zeros(size(F.x)) ;
    F.ab=zeros(size(F.x)) ;
    F.dasdh=zeros(size(F.x)) ;
    F.dabdh=zeros(size(F.x)) ;
end




dabdh=F.ab;
ab=F.dabdh;


I=F.b>dMin ; ab(I)=abMin; dabdh(I)=abMin;        % above dMin

I= F.b< dMin & F.b >= dMax ; ab(I)=abMax*(F.b(I)-dMin)/(dMax-dMin); 
dabdh(I)=(-F.rho(I)/F.rhow)/(abMax/(dMax-dMin));

I=F.b<dMax ; ab(I)=abMax; dabdh(I)=abMax;        % below dMax




end