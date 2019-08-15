
function [L,cuvh,luvh]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs,l)

MLC=BCs2MLC(CtrlVar,MUA,BCs);
Luv=MLC.ubvbL;
cuv=MLC.ubvbRhs;
Lh=MLC.hL;
ch=MLC.hRhs;


% Luv  : #uv constrains x 2MUA.Nnodes
% Lh  : #h constrains x MUA.Nnodes
% L=[Luv 0]
%   [0  Lh]

[nu,~]=size(Luv) ;
[nh,~]=size(Lh) ;


if CtrlVar.LinFEbasis
    
    % L M L' L
    %
    % L -> L M
    % c -> L M L' c
    if ~isfield(MUA,'M')
        MUA.M=MassMatrix2D1dof(MUA);
    end
    
    Mblock=MassMatrixBlockDiagonal2D(MUA);
    
    if numel(Luv)>0
        Luv=Luv*Mblock ;
    end
    
    if numel(Lh)>0
        Lh=Lh*MUA.M  ;
    end
    
    if numel(cuv)>0
        cuv=(Luv*Mblock*Luv')*cuv ;
    end
    
    if numel(ch)>0
        ch=(Lh*MUA.M*Lh')*ch ;
    end
    
end


if isempty(Lh) && ~isempty(Luv)
    L=[Luv sparse(nu,MUA.Nnodes)] ;
    cuvh=cuv ;
    luvh=l.ubvb;
elseif ~isempty(Lh) && isempty(Luv)
    L=[sparse(nh,2*MUA.Nnodes) Lh] ;
    cuvh=ch ;
    luvh=l.h;
elseif ~isempty(Lh) && ~isempty(Luv)
    
    L=[ Luv sparse(nu,MUA.Nnodes) ; sparse(nh,2*MUA.Nnodes) Lh];
    cuvh=[cuv;ch];
    luvh=[l.ubvb;l.h];
else
    L=[] ; cuvh=[] ; luvh=[];
end


end

