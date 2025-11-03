function [ufixedvalue,vfixedvalue] =GetPIGfixedVelocities(DirichletNodes,coordinates)
	
	CurDir=pwd;
	goto_home_directory;
	cd('PIG')
	load PIG96VelocityHilmar X Y Speed Vx Vy
%	load BamberVelDataGridded X Y Vx Vy Speed 	;   X=1000*X ; Y=1000*Y ; % go from km to m
	cd(CurDir)
	
	%%
	
	% somewhat questionable, but will have to do for the time being
	Vx(isnan(Vx))=0;
	Vy(isnan(Vy))=0;
	
	ufixedvalue=interp2(X,Y,Vx,coordinates(DirichletNodes,1),coordinates(DirichletNodes,2),'linear');
	vfixedvalue=interp2(X,Y,Vy,coordinates(DirichletNodes,1),coordinates(DirichletNodes,2),'linear');
	
	
	
end

