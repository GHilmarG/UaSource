function [AGlen,C,CAGlen,CC,lxfixednode,lyfixednode,lxfixedvalue,lyfixedvalue,lxtiedA,lxtiedB,lytiedA,lytiedB...
		sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,xMeas,yMeas,DTmeas]=...
		DefineInverseModellingVariables(AGlen,C,s,b,B,ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,CtrlVar)
	
	
	if strcmp(Experiment,'inverse1') || strcmp(Experiment,'forwardC')|| strcmp(Experiment,'forwardA')
		
		CError=zeros(Nnodes,1)+C0; CC=sparse(1:Nnodes,1:Nnodes,1./CError.^2,Nnodes,Nnodes);   % I=(B u-d)' M (B u -d )
		AGlenError=zeros(Nnodes,1)+mean(AGlen) ; CAGlen=sparse(1:Nnodes,1:Nnodes,1./AGlenError.^2,Nnodes,Nnodes);   % I=(B u-d)' M (B u -d )
		
		
		AGlen=AGlen+zeros(Nnodes,1)  ; C=zeros(Nnodes,1)+C0;
		
		
		fprintf(' loading synthetic data for inverse step and applying Dirchlet BCs for inverse step \n')
		load(CtrlVar.SurfaceDataFile,'sMeas','uMeas','vMeas','wMeas','bMeas','BMeas','xMeas','yMeas')
		CtrlVar.SurfaceDataFileLoaded=1;
		DTmeas = DelaunayTri(xMeas,yMeas); 
		
		
		switch lower(CtrlVar.AdjointBC)
			case 'periodic'
				lxtiedA=utiedA ; lxtiedB=utiedB;
				lytiedA=vtiedA ; lytiedB=vtiedB;
				lxfixednode=[];  lxfixedvalue=[];
				lyfixednode=[];  lyfixedvalue=[];
			case 'dirichletzero'
				lxfixednode=Boundary.Nodes;   lyfixednode=Boundary.Nodes;
				lxfixedvalue=lxfixednode*0 ;  lyfixedvalue=lyfixednode*0 ;
				lxtiedA=[]; lxtiedB=[]; lytiedA=[] ; lytiedB=[];
			otherwise
				error(' case not recognised ')
		end
		
		
	elseif strcmp(Experiment,'PIG')
		
		% for inverse step:
		if CtrlVar.doInverseStep==1
			fprintf(' Loading measured velocity data for PIG \n')
			BMeas=B ; bMeas=b ; sMeas=s ;
			CurDir=pwd; goto_home_directory ; cd('PIG') ;  load PIG96VelocityHilmar X Y Speed Vx Vy ;cd(CurDir)
			[uMeas,vMeas] = PIGfVelocity(X,Y,Vx,Vy,coordinates,1:Nnodes);
			wMeas=uMeas*0;
			xMeas=coordinates(:,1); yMeas=coordinates(:,2);
			DTmeas = DelaunayTri(xMeas,yMeas);
		end
		clear fdata
		
		
	end
	
	
end
