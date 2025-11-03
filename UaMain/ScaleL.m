function [L,c,isLLunity,Scale]=ScaleL(CtrlVar,L,c,test)

%%
%
% Scales L so that L*L'=1
%
% Provided L has only one non-zero element in each of its columns, this will
% always work.
%
% If L is a MLC matrix where there is only one constraint for each constrained node,
% L will only have one non-zero element along each column. Hence, for a properly
% constructed multi-linear constraints, this scaling will ensure that
% L*L=ones(p,p), where p is the number of lines in L, i.e. the number of
% constraints.
%
%%

isLLunity=[];

if isempty(L)
    return
end

if nargin<4
    test=false ;
end

p=size(L,1);
Lscale=sqrt(sum(abs(L).^2,2));
Scale=spdiags(1./Lscale,0,p,p);

L=Scale*L ;
c=Scale*c ;


% isLLunity=norm(L*L'-eye(p,p))< 100*eps ;  % slow
if test
    isLLunity=full(sum(sum(L*L'-speye(p,p)))) < 100*eps ;
end


end