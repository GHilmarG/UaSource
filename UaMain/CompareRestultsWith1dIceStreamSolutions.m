function CompareRestultsWith1dIceStreamSolutions(Experiment,coordinates,connectivity,u,v,s,b,S,B,time,dt,AGlen,C,n,m,rho,rhow,g,alpha,nip,DTxy,TRIxy,DTint,TRIint,CtrlVar)
    
    
    
    
    x=coordinates(:,1);
    
    GF = GL2d(B,S,s-b,rhow,rho,connectivity,CtrlVar);  % gf is only needed for plotting purposes and global remeshing
    
    
    [GLpos,indGL]=max(x.*GF.node);
    
    
    figure(1) ; hold on 
    
    [Nele,nod]=size(connectivity);
    
    switch nod
        case 3
            col='b';
        case 6
            col='k';
        case 10
            col='g';
    end
    
    figure(1001)
    hold on
    plot(x/1000,u,'x','color',col)
    hold on ; 
    %plot([ GLpos/1000 GLpos/1000],[min(u) max(u)])
    xlabel('km') ; title(' u ')
    
    alpha=0.001; dGL=S(indGL)-b(indGL);
    [uAna,xAna,GLposAna,dudxGLAna,dl]=Simple1dIceStreamSolution(x,s,b,S,B,alpha,g,rho,rhow,AGlen,C);
    
    
    
    fprintf('\t  GLposAna=%-g km\t dl=%-g m\t dudxGLAna=%-g \n ',GLposAna/1000,dl,dudxGLAna)
  %  fprintf('\t  GLpos=%-g km   \t dGL=%-g m \n ',GLpos/1000,dGL);
  
    
    plot(xAna/1000,uAna,'r')
    

end
