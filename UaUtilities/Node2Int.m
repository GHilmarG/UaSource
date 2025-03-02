

function  [xInt,yInt,varargout]=Node2Int(CtrlVar,MUA,varargin)

%%
%
% Calculates integration point values, and the (x,y) coordinates of the integration points
%
% Example:
%
% 
%   [xInt,yInt,h0int,B0int,S0int,rho0int]=Node2Int(CtrlVar,MUA,F0.h,F0.B,F0.S,F0.rho);
%   UaPlots(CtrlVar,MUA,F1,h0int)
%
%
%
%%


nVarargsIn = length(varargin);
nVar=nVarargsIn;

fnod=cell(nVar,1);

ndim=2;

coox=reshape(MUA.coordinates(MUA.connectivity,1),MUA.Nele,MUA.nod);
cooy=reshape(MUA.coordinates(MUA.connectivity,2),MUA.Nele,MUA.nod);

for I=1:nVarargsIn
    fnod{I}=reshape(varargin{I}(MUA.connectivity,1),MUA.Nele,MUA.nod);
end

xInt=zeros(MUA.Nele,MUA.nip) ;
yInt=zeros(MUA.Nele,MUA.nip);

varargout=cell(nVar,1);

for I=1:nVar
    varargout{I}=zeros(MUA.Nele,MUA.nip) ;
end

% vector over all elements for each integration point

for Iint=1:MUA.nip


    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points

    xInt(:,Iint)=coox*fun;
    yInt(:,Iint)=cooy*fun;
    for I=1:nVar
        varargout{I}(:,Iint)=fnod{I}*fun;
    end

end

end

