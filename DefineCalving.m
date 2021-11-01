function [UserVar,LSF,c]=DefineCalving(UserVar,CtrlVar,MUA,F,BCs)
    

%%
%
%   [UserVar,LSF,CalvingRate]=DefineCalving(UserVar,CtrlVar,MUA,F,BCs)
%
% Define calving the Level-Set Field (LSF) and the Calving Rate Field (c)
%
%
% NOTE: Currently the only supported calving option is based on directly prescribing the LSF at each time step
%
% Think of the LSF as a ice/ocean mask where positive number indicates ice and negative number ocean.
%
% LSF is a nodal variable.
%
% If you wanted, for example to get rid of all floating ice for x>500e3, do:
%
%
%
%
%         F.GF=IceSheetIceShelves(CtrlVar,MUA,F.GF);
%         OceanNodes=MUA.coordinates(:,1)>500e3 & F.GF.NodesDownstreamOfGroundingLines;
%         LSF=zeros(MUA.Nnodes,1)+ 1 ;
%         LSF(OceanNodes)=-1;
%
%
%     
%%

LSF=[];
c=[];

end
