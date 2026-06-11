function [R,K,F,T]=MatrixAssemblySSHEETtransient2HD(CtrlVar,MUA,AGlen,n,C,m,rho,g,h0,b0,h1,b1,a0,a1,dt)

narginchk(15,15)

if nargout==1  
    OnlyR=1 ; 
else
    OnlyR=0 ;
end

if CtrlVar.InfoLevelCPU>=10 ;   tCPU=tic; end

% NR:  h^{i+1} = \Delta h + h_1  and h^i=h_0
%      I solve for \Delta h.
%      If the thing converges then \Delta h -> 0 as i->\infty

% h0=s0-b0;
% h1=s1-b1;

s0=h0+b0;
s1=h1+b1;

neq=MUA.Nnodes;



h0nod=reshape(h0(MUA.connectivity,1),MUA.Nele,MUA.nod);
h1nod=reshape(h1(MUA.connectivity,1),MUA.Nele,MUA.nod);
s0nod=reshape(s0(MUA.connectivity,1),MUA.Nele,MUA.nod);
s1nod=reshape(s1(MUA.connectivity,1),MUA.Nele,MUA.nod);
a0nod=reshape(a0(MUA.connectivity,1),MUA.Nele,MUA.nod);
a1nod=reshape(a1(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

AGlen=reshape(AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
n=reshape(n(MUA.connectivity,1),MUA.Nele,MUA.nod);

C=reshape(C(MUA.connectivity,1),MUA.Nele,MUA.nod);
m=reshape(m(MUA.connectivity,1),MUA.Nele,MUA.nod);


% can I not just use h0 for h0MUA.nod, ie reuse variable?


K=sparse(neq,neq);

if ~OnlyR
    d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
else
    d1d1=[];
end

t1=zeros(MUA.Nele,MUA.nod); f1=zeros(MUA.Nele,MUA.nod);


% vector over all elements for each integartion point

if CtrlVar.Parallel.hAssembly.parfor.isOn
    poolobj = gcp;
    Mworkers=poolobj.NumWorkers;
else
    Mworkers=1;
end

%parfor (Iint=1:MUA.nip,Mworkers)
for Iint=1:MUA.nip
    
    [tt,ff,dddd]=SSHEETintAssembly(Iint,CtrlVar,MUA,AGlen,n,C,m,rhonod,g,s0nod,h0nod,s1nod,h1nod,a0nod,a1nod,dt,OnlyR);
    t1=t1+tt;
    f1=f1+ff;
    d1d1=d1d1+dddd;
end


T=sparseUA(neq,1);
for I=1:MUA.nod
    T=T+sparseUA(MUA.connectivity(:,I),ones(MUA.Nele,1),t1(:,I),neq,1);
end

F=sparseUA(neq,1);
for I=1:MUA.nod
    F=F+sparseUA(MUA.connectivity(:,I),ones(MUA.Nele,1),f1(:,I),neq,1);
end

R=T-F  ;

if ~OnlyR
    for I=1:MUA.nod
        for J=1:MUA.nod
            K=K+sparseUA(MUA.connectivity(:,I),MUA.connectivity(:,J),d1d1(:,I,J),neq,neq);
        end
    end
else
    K=[];
end

%K=(K+K')/2;

if CtrlVar.InfoLevelCPU>=10 ;  fprintf(CtrlVar.fidlog,' MatrixAssemblySSHEETtransient2HD in %-g sec. \n',toc(tCPU)) ; end



end