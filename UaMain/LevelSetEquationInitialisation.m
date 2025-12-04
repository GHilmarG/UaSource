function [UserVar,RunInfo,LSF,l,LSFqx,LSFqy,BCs]=LevelSetEquationInitialisation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)

nargoutchk(7,7)

LSFqx=[]; LSFqy=[] ; 

fprintf("LevelSetEquationInitialisation: Initialising the level set. \n ")

switch lower(CtrlVar.LevelSetInitialisationMethod)

    

    case {"-geometric-","geometric","-geo-","geo"}

        Value=0 ;  [Xc,Yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F0.LSF,Value,subdivide=true) ;
        isPlot=false;
        CFPD=sqrt(2*min(MUA.EleAreas))/CtrlVar.LevelSetGeometricInitialisationDistanceFactor;
        [~,~,LSF]=...
            CalvingFrontLevelSetGeometricalInitialisation(CtrlVar,MUA,Xc,Yc,F0.LSF,...
            method="InputPoints",...
            ResampleCalvingFront=true,...
            CalvingFrontPointDistance=CFPD,...
            plot=isPlot) ;
        %%

    otherwise

        % CtrlVar.LevelSetReinitializePDist=false ;
        [UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationAnalyticalInitialisation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
        F0.LSF=LSF ; F1.LSF=LSF ;
end






end





