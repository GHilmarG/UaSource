
mem=0.25;

matlabpool(8)
results8=paralleldemo_backslash_benchGHG(mem);
matlabpool close


matlabpool(4)
results4=paralleldemo_backslash_benchGHG(mem);
matlabpool close


matlabpool(2)
results2=paralleldemo_backslash_benchGHG(mem);
matlabpool close

matlabpool(1)
results1=paralleldemo_backslash_benchGHG(mem);
matlabpool close

%%
figure
plot(results8.matSize,results8.gflops,'-+b')
hold on
plot(results4.matSize,results4.gflops,'-+r')
plot(results2.matSize,results2.gflops,'-+g')
plot(results1.matSize,results1.gflops,'-+c')

legend('8','4','2','1')
xlabel('Matrix size')
ylabel('Gflops')
