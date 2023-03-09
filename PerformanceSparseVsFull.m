
%%




% density=0.001;
% iRepeat=10;
% 
% 
% fprintf('\n \n N=%i \n ', N)
% fprintf(' Creating arrays.\n')
% 
% Asparse=sprand(N,N,density) ;
% Asparse=Asparse+sparse(1:N,1:N,1) ;
% Afull=full(Asparse);
% x=ones(N,1) ;


% If not available locally, download from
% https://livenorthumbriaac-my.sharepoint.com/:f:/g/personal/hilmar_gudmundsson_northumbria_ac_uk/EqKCbhYlStJKlSCUBnq-n8MBHLY-oJbdRovjZWMHyqJU9Q?e=IRG2Iz
load("ABmatrixExamplesA243858B2568.mat","A","B","f","g")


% 2022-09-22: C17777347 with Quadro M4000
% N=40000   full=NaN  	 sparse=0.754783 	 gpu=3.402497 	 full/sparse=NaN gpuSparse/cpuSparse=4.507915 

doFull=false; 
N=50000; iRepeat=4; 
Asparse=A(1:N,1:N);
Afull=full(Asparse);
AGPUsparse=gpuArray(Asparse) ;

f=f(1:N);

if doFull
    fprintf("CPU full \n")
    CPUfull=tic ;
    for I=1:iRepeat
        yfull=Afull\f;
    end
    CPUfull=toc(CPUfull)/iRepeat;
else
    CPUfull=nan;
end


fprintf("CPU sparse \n")
CPUsparse=tic ;
for I=1:iRepeat
    yfull=Asparse\f;
end
CPUsparse=toc(CPUsparse)/iRepeat;

fprintf("GPU sparse \n")
gd = gpuDevice();
GPUsparse=tic ;
for I=1:iRepeat
    yfull=AGPUsparse\f;
end
wait(gd); 
GPUsparse=toc(GPUsparse)/iRepeat;


fprintf("full=%f  \t sparse=%f \t gpu=%f \t full/sparse=%f gpuSparse/cpuSparse=%f \n",CPUfull,CPUsparse,GPUsparse,CPUfull/CPUsparse,GPUsparse/CPUsparse)






