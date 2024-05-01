function Degree=QuadratureRuleDegree(CtrlVar)



if isfield(CtrlVar,"QuadratureRuleDegree")  && ~isempty(CtrlVar.QuadratureRuleDegree) && ~isnan(CtrlVar.QuadratureRuleDegree)
    
    Degree=CtrlVar.QuadratureRuleDegree;

else

    switch CtrlVar.TriNodes

        case 3

            Degree=4;
            Degree=25;


        case 6

            Degree=8;
           % Degree=12;
           % Degree=20;
           % Degree=25;

        case 10

            Degree=8 ;
            % Degree=25;

        otherwise

            error('Ua:CaseNotFound','Which case?')

    end

end


end