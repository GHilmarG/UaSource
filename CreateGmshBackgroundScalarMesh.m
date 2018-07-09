function io = CreateGmshBackgroundScalarMesh(arg1,arg2,TargetSize,FileName)
    % creates a Gmsh background mesh file that can be used to define
    % target element sizes at location of nodal points (xy)
    % The output is written to a file named
    %
    % arg1 and arg2 can be either:
    % xTRI and yTRI arrays of size N x 3 defining the xy coordinates of 3 node triangles
    % or
    % coordinates (nNodes x 2 ) and connectivity (nEle x 3) of 3 node triangle mesh
    %
    % The output file only contains the coordinats of the corners of the (linear) triangles
    % the connectivity as such does not appear in that file, and the output file is independent of
    % nodal numbering


    io=0;
    [n1,m1]=size(arg1) ;
    [n2,m2]=size(arg2) ;

   
   

    if m1==3 && m2==3
        xTRI=arg1 ; yTRI=arg2 ;
        
    elseif m1==2 && m2==3
        
        x=arg1(:,1);  y=arg1(:,2); connectivity=arg2;
        xTRI=x(connectivity) ; yTRI=y(connectivity) ; TargetSize=TargetSize(connectivity);
        
    else
        
        error('inputs not consistent')
        io=1;
    end
    
    z=zeros(length(xTRI),1);
    
    xyzTRI=[xTRI(:,1) yTRI(:,1) z xTRI(:,2) yTRI(:,2) z xTRI(:,3) yTRI(:,3) z];
    
    Temp=[xyzTRI TargetSize];
%     figure
%     plot3(Temp(:,1),Temp(:,2),Temp(:,10),'o') ; hold on
%     plot3(Temp(:,4),Temp(:,5),Temp(:,10),'o')
%     plot3(Temp(:,7),Temp(:,8),Temp(:,10),'o') ; hold off
%     figure
%     
    fileID = fopen(FileName,'w');
    
    Temp=Temp';

    fprintf(fileID,'View "background mesh" {\n');
    fprintf(fileID,'ST(%f,%f,%f,%f,%f,%f,%f,%f,%f){%f,%f,%f};\n',Temp);
    fprintf(fileID,'};');
    fclose(fileID);
    
   
    
    
end

