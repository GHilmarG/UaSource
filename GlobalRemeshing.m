
function  [JobVar,coordinates,connectivity,...
        ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,hfixednode,hfixedvalue,htiedA,htiedB,...
        Luv,Luvrhs,lambdauv,Lh,Lhrhs,lambdah,...
        Boundary,u,v,h,s,b,S,B,AGlen,n,C,m,rho,alpha,DTxy,TRIxy,DTint,TRIint,Xint,Yint,dhdt,dudt,dvdt,dhdtm1,dudtm1,dvdtm1,...
        GF,xGLmesh,yGLmesh]=...
        GlobalRemeshing(JobVar,Experiment,MeshBoundaryCoordinates,Boundary,...
        s,b,S,B,h,u,v,wSurf,coordinates,connectivity,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,time,Itime,CtrlVar,...
        dhdt,dudt,dvdt,dhdtm1,dudtm1,dvdtm1,DTxy,TRIxy)
    
   xGLmesh=[] ; yGLmesh=[];
 
    
    
    % remesh domain and interpolate values onto the new grid
    
    %%
    
    
    % save TestSave;  error('dsfa')
    
    switch lower(CtrlVar.MeshRefinementMethod)
        
        case 'implicit'
            
            [JobVar,coordinates,connectivity]=DesiredEleSizesBasedOnImplicitErrorEstimate(JobVar,Experiment,MeshBoundaryCoordinates,Boundary,...
                s,b,S,B,h,u,v,coordinates,connectivity,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,DTxy,TRIxy,CtrlVar);
  
        case 'explicit'
            
            [JobVar,coordinates,connectivity,xGLmesh,yGLmesh]=...
                DesiredEleSizesBasedOnExplicitErrorEstimate(JobVar,Experiment,MeshBoundaryCoordinates,...
                S,B,h,s,b,u,v,wSurf,dhdt,coordinates,connectivity,nip,AGlen,C,n,rho,rhow,DTxy,TRIxy,CtrlVar);
            
        otherwise
            
            error(' unknown case  ')
    end

  
    
end
