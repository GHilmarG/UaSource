

function [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

%%
%  User input m-file to define A and n in the Glenn-Steinemann flow law
%
%   [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)
%
% Usually A is defined on the nodes, but sometimes in an inverse run A might be
% defined as an element variable. The user makes this decision by setting
% CtrlVar.AGlenisElementBased to true or false in Ua2D_InitialUserInput
%%
warning('Ua:DefaultDefine','Using default DefineAGlenDistribution')


n=3 ;

if CtrlVar.AGlenisElementBased
    AGlen=3.0e-9+zeros(MUA.Nele,1) ;  % kPa year about -20 degrees Celcius
else
    AGlen=3.0e-9+zeros(MUA.Nnodes,1) ; % kPa year about -20 degrees Celcius
end



if CtrlVar.doDiagnostic
    
    switch lower(UserVar.RunType)
        
        case 'icestream'
            
        case 'iceshelf'
            x=MUA.coordinates(:,1) ;
            y=MUA.coordinates(:,2);
            
            if CtrlVar.CisElementBased
                x=mean(reshape(x(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
                y=mean(reshape(y(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
            end
            
            sx=10e3 ; sy=10e3;
            AGlen=AGlen.*(1+100*exp(-(x.*x/sx^2+y.*y./sy^2)));

            
    end
end





end

