
%%
N=100doc  ; density=0.1 ; 

A=sprandsym(N,density) ; x=ones(N,1) ; 
%A=rand(N,N) ; x=ones(N,1) ; 

%Agpu=gpuArray(single(A));  xgpu=gpuArray(single(x)) ;
Agpu=gpuArray(A);  xgpu=gpuArray(x) ;

fprintf('\n\n')
tic 
y=A\x;
toc

tic
ygpu=Agpu\xgpu;
toc

