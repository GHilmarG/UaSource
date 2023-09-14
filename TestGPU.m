
function TestGPU
%%
%
% conclusion on 19 April,2019. Same conclusions on 14 Sep, 2023.
% GPU linsolve using sparse matrices about 3 to 4 times slower than CPU
% GPU linsolve using full matrices a bit faster than CPU
%
% 
%%

iExperiment=0;
density=0.05 ;
timings=zeros(10,7)+NaN;

nRepeat=2;

for N=[100 1000 5000 10000] % 500 1000 2000 3000 20000]
    
    
    iExperiment=iExperiment+1;
    
    fprintf('\n \n N=%i \n ', N)
    fprintf(' Creating arrays.\n')
    
    Asparse=sprand(N,N,density) ;
    Asparse=Asparse+sparse(1:N,1:N,1) ; 
    Afull=full(Asparse);
    x=ones(N,1) ;
    Adistfull=distributed(Afull) ;
    Adistsparse=distributed(Asparse) ;
    Agpufull=gpuArray(Afull);
    xgpu=gpuArray(x) ;
    Agpusparse=gpuArray(Asparse) ;
    xdist=distributed(x);
    
    fprintf('A CPU distributed full. \n')
    for k=1:nRepeat
        CPUdistfull=tic ;
        y=Adistfull\xdist;
        CPUdistfull=toc(CPUdistfull);
    end
    
    CPUdistsparse=NaN;
    %% matlab 2019a does not work with sprandsym, only sprand
    fprintf('A CPU distributed sparse. \n')
    for k=1:nRepeat
        CPUdistsparse=tic ;
        y=Adistsparse\xdist;
        CPUdistsparse=toc(CPUdistsparse);
    end
    %%
    
    fprintf('A CPU full  \n')
    for k=1:nRepeat
        CPUfull=tic ;
        y=Afull\x;
        CPUfull=toc(CPUfull);
    end
    
    fprintf('A CPU sparse \n')
    for k=1:nRepeat
        CPUsparse=tic ;
        y=Asparse\x;
        CPUsparse=toc(CPUsparse);
    end
    
    fprintf('A GPU \n')
    for k=1:nRepeat
        
        GPUfull=tic;
        ygpu=Agpufull\xgpu;
        GPUfull=toc(GPUfull);
    end

    fprintf('A GPU sparse \n')
    for k=1:nRepeat
        GPUsparse=tic ;
        y=Agpusparse\xgpu;
        GPUsparse=toc(GPUsparse);
    end

    %     Agpu=single(gpuArray(A));
    %     xgpu=single(gpuArray(x)) ;
    %     GPUsingle=tic;
    %     ygpu=Agpu\xgpu;
    %     GPUsingle=toc(GPUsingle);

    timings(iExperiment,1)= N ;  timings(iExperiment,2)= CPUfull ;
    timings(iExperiment,3)= GPUfull ;
    timings(iExperiment,4)= CPUdistfull ;
    timings(iExperiment,5)= CPUsparse ;
    timings(iExperiment,6)= CPUdistsparse ;
    timings(iExperiment,7)= GPUsparse ;
    
end

figure
plot(timings(:,1),timings(:,2),'o-r')
hold on
plot(timings(:,1),timings(:,3),'x-b')
plot(timings(:,1),timings(:,4),'+-g')
plot(timings(:,1),timings(:,5),'*-c')
plot(timings(:,1),timings(:,6),'^-m')
plot(timings(:,1),timings(:,7),'o-k')
legend('CPU full','GPU double full','CPU distributed full','CPU sparse','CPU dist sparse','GPU sparse')
xlabel("Problem size N")
ylabel("time (sec)")

end