function [uError,vError,wError]=DefineDataErrors(coordinates,Xint,Yint,CtrlVar)
    
    
    Nnodes=length(coordinates);
     
    uError=zeros(Nnodes,1)+1; vError=zeros(Nnodes,1)+1;  nInt=numel(Xint) ;  wError=zeros(nInt,1)+1;
    
end
