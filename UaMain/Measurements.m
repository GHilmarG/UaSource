classdef Measurements

    properties


        us=[];
        vs=[];
        dhdt=[] ;


        usCov=[];
        vsCov=[];
        dhdtCov=[];

        s=[] ; sCov=[];

        Bobs=[];   % Direct observations of B
        Bx=[];     % x location of direct observations
        By=[];     % y location of direct observations
        BErr=[];   % Error, actually this is sigma
        BO=[];      % Node2Data mapping matrix, this is a property of both the data and the mesh!
        BInside=[]; % logical array indicating if measurement i is inside or outside the mesh
        BEleID=[];  % element number of the i-th measurement. 


        as=[] ;
        ab=[] ; 


    end

end