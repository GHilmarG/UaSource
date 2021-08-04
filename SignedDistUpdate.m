function [LSF,UserVar,RunInfo]=SignedDistUpdate(UserVar,RunInfo,CtrlVar,MUA,LSF,xc,yc)


narginchk(7,7)


if numel(xc)>0
    
    Dist=pdist2([xc(:) yc(:)],MUA.coordinates,'euclidean','Smallest',1) ;
    Dist=Dist(:) ;

    PM=sign(LSF) ;
    LSF=PM.*Dist;
    
else
    
    fprintf('SingedDistUpdate:No nodes on input.\n')
    
end



end