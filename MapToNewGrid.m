function [s,b,h,S,B,rho,AGlen,n,C,m,GF,DTxy2,TRIxy2,DTint2,TRIint2,Xint2,Yint2,xint2,yint2,Iint2,varargout]=...
        MapToNewGrid(CtrlVar,MUA,DTxyInput,DTxy1,hInput,Itime,time,varargin)

    
    
    Nin=7 ;Nout=20;
    
    if (nargin-Nin) ~= (nargout-Nout)
        error(' incorrect combination of input and output variables ')
    end
    
    x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
    
    if CtrlVar.doDiagnostic
        [UserVar,s,b,S,B,alpha]=DefineGeometry(CtrlVar.Experiment,MUA.coordinates,CtrlVar,'sbSB');
        h=s-b;
    else
        [~,~,S,B,alpha]=DefineGeometry(CtrlVar.Experiment,MUA.coordinates,CtrlVar,'SB');
        h=Grid1toGrid2(DTxyInput,hInput,x,y);
    end

    
    %[UserVar,rho,rhow,g]=DefineDensities(Experiment,coordinates,connectivity,s,b,h,S,B,Itime,time,CtrlVar)
    [UserVar,rho,rhow,g]=DefineDensities(CtrlVar.Experiment,MUA.coordinates,MUA.connectivity,[],[],h,S,B,Itime,time,CtrlVar); 
    [b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);

    [UserVar,C,m]=DefineSlipperyDistribution(CtrlVar.Experiment,MUA.coordinates,MUA.connectivity,s,b,h,S,B,rho,rhow,Itime,time,CtrlVar);
    [UserVar,AGlen,n]=DefineAGlenDistribution(CtrlVar.Experiment,MUA.coordinates,MUA.connectivity,s,b,h,S,B,rho,rhow,Itime,time,CtrlVar);
    
    GF = GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);  
    [DTxy2,TRIxy2,DTint2,TRIint2,Xint2,Yint2,xint2,yint2,Iint2]=TriangulationNodesIntegrationPoints(MUA);
    
    
    
    
    nIn=nargin-Nin;
    if nargout>Nout
        varargout=cell(nargout-Nout,1);
        for I=1:nIn
            varargout{I}=Grid1toGrid2(DTxy1,varargin{I},x,y);
        end
    end
    
    
    
    
    
end