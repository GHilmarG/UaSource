

%%
ndim=2;

nipv=  [ 1,3,4,6,7, 9,12,13,16,19,28,37];
degree=[ 1 2 3 4 5 6  6  7  8  9 11 13];
figure ;
I=0;
for nip=nipv
    [s,wt] = sample('triangle',nip,ndim);
    I=I+1;
    subplot(3,4,I)
    plot([0 0 1 0],[0 1 0 0],'LineWidth',2)
    hold on
    plot(s(:,1),s(:,2),'o','MarkerFaceColor','r')
    axis equal
    title(sprintf(' nip=%i, degree=%i ',nip,degree(I)))
    xlabel('x'); ylabel('y'); 
end