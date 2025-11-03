function [UserVar,R,K,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistent(UserVar,CtrlVar,MUA,F0,F1)

% [UserVar,R,K,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistent(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1)


narginchk(5,5)
nargoutchk(2,6)


switch CtrlVar.LevelSetMethodEquationForm

    case "scalar"

        [UserVar,R,K]=LevelSetEquationAssemblyNR2consistentScalar(UserVar,CtrlVar,MUA,F0,F1) ;

        % [UserVar,R,K]=LevelSetEquationAssemblyNR2consistentScalar(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1);
        Qx=[] ; Qy=[] ; Rv=[];

    case "vector"

        error("NotUpdatedToReflectLaterChanges")
        [UserVar,R,K,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistentVector(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1);

    otherwise

        error("LevelSetEquationAssemblyNR2consistent:CaseNotFound","Case not found")

end



end