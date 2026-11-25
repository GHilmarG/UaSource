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
        BLowerBound=[];   % (20 Nov 2026) optional node-wise lower bound on B for the inversion (NaN: none), see F2p.m
        BUpperBound=[];   % (20 Nov 2026) optional node-wise upper bound on B for the inversion (NaN: none), see F2p.m
        BO=[];      % Node2Data mapping matrix, this is a property of both the data and the mesh!
        BInside=[]; % logical array indicating if measurement i is inside or outside the mesh
        BEleID=[];  % element number of the i-th measurement. 

        % Point measurements of surface velocities, used if CtrlVar.Inverse.DataMisfit.uv.Method="-points-" (see CreatePointMeasurements.m)
        uvx=[];     % x locations of the velocity measurements
        uvy=[];     % y locations of the velocity measurements
        usData=[];  % measured u at (uvx,uvy)
        vsData=[];  % measured v at (uvx,uvy)
        usErr=[];   % errors (standard deviations) of usData
        vsErr=[];   % errors (standard deviations) of vsData
        uvO=[];     % Node2Data mapping matrix for (uvx,uvy), this is a property of both the data and the mesh!
        uvCorr=[];  % (23 Nov 2026) spatially correlated velocity errors: [] (default, errors usErr, vsErr uncorrelated) or a structure with
                    % fields sigma and rho (Matern, alpha=2, on the mesh); usErr, vsErr are then the uncorrelated parts (nugget), see
                    % CorrelatedPointDataMisfit.m

        % Measurements of dh/dt used if CtrlVar.Inverse.DataMisfit.dhdt.Method="-points-" or "-blocks-"
        % (see CreateDhdtPointMeasurements.m and CreateDhdtBlockMeasurements.m)
        dhdtx=[];          % x locations of the measurements ("-points-"), or of the block centres ("-blocks-")
        dhdty=[];          % y locations
        dhdtData=[];       % measured dh/dt ("-points-"), or block averages of measured dh/dt ("-blocks-")
        dhdtErr=[];        % errors (standard deviations) of dhdtData
        dhdtP=[];          % ("-points-") form functions and derivatives at the measurement locations (ShapeFunctionsAtPoints.m), depends on the mesh!
        dhdtA=[];          % ("-blocks-") block-averaging operator acting on values at the integration points, depends on the mesh!
        dhdtBlockSize=[];  % ("-blocks-") block size


        as=[] ;
        ab=[] ; 


    end

end