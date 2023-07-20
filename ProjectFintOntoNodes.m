function varargout=ProjectFintOntoNodes(MUA,varargin)

%%
%   varargout=ProjectFintOntoNodes(MUA,varargin)
%
% Projects element variables defined at integration points onto nodes.
% Does this by solving:
%
%          min_{Fnod_p} || Fnod_p N_P - Fint||_{L2}
%
%  ->  <Fnod_p N_p - Fint , N_q >_{L2}=0
%  ->  <N_p,N_q> Fnod_p = <Fint,N_q>
%        M Fnod = <Fint,n_q>
%
% No limit on number of input/output fields.
%
% Example:
%
%  [exxNod,eyyNod,exyNod,eNod]=ProjectFintOntoNodes(MUA,exx,eyy,exy,e)
%
% where exx, eyy, exy, e, are defined at integration points gives
% corresponding fields defined at nodes.
%
% Note: This projection will not always preserve positivity! 
%       Even if element quantity is positive everywere, it is nevertheless 
%       possible that the projection on the nodes can be (slightly) negative.
%       


nVarargs = length(varargin);
varargout = cell(nVarargs);


% check input dimentions
for I=1:nVarargs
    [N,M]=size(varargin{I});
    if N~=MUA.Nele || M~=MUA.nip
        fprintf('Incorrect dimensions: Must be an element variable defined at all elements and all integration points\n')
        return
    end
end

% create mass matrix

b=zeros(MUA.Nnodes,nVarargs);

% factorize
% [L,~,P]=chol(A,'lower');
for I=1:nVarargs
    b(:,I)=InnerProduct_FormFunctions_with_EleIntegrationPointVariable(MUA,varargin{I});
    %      varargout{I}= P*(L' \(L \(P'*b(:,I))));
end

%A=MassMatrix2D1dof(MUA);
%sol=A\b;

if ~isfield(MUA,'M') || isempty(MUA.M)
    MUA.M=MassMatrix2D1dof(MUA);
end

if isfield(MUA,"dM") && isa(MUA.dM,"decomposition") && ~isempty(MUA.dM)
    try
        sol=MUA.dM\b;
    catch
        sol=MUA.M\b;
    end
else
    sol=MUA.M\b;
end

for I=1:nVarargs
    varargout{I}=sol(:,I);
end

end

