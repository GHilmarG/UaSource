function Degree=QuadratureRuleDegree(CtrlVar)



if isfield(CtrlVar,"QuadratureRuleDegree")  && ~isempty(CtrlVar.QuadratureRuleDegree) && ~isnan(CtrlVar.QuadratureRuleDegree)
    
    Degree=CtrlVar.QuadratureRuleDegree;

else

    switch CtrlVar.TriNodes

        case 3

            Degree=4;


        case 6

            Degree=8;
           % Degree=12;
           % Degree=20;

        case 10

            Degree=8 ;

        otherwise

            error('Ua:CaseNotFound','Which case?')

    end

end


end