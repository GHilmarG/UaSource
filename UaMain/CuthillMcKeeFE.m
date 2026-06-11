function [coo,connectivity] = CuthillMcKeeFE(coo,connectivity,M)
    
    %[M] = connectivity2adjacency(connectivity);
    p=symrcm(M);  disp(' Reverse Cuthill-McKee')
    %p=symamd(M);
    
    p2=p*0; p2(p)=1:length(p2);
    coo(:,2)=coo(p,2);
    coo(:,1)=coo(p,1);
    connectivity=p2(connectivity);
    
    
end

