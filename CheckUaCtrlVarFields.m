
function CheckUaCtrlVarFields(CtrlVar,Field)

if ~ismember(string(CtrlVar.(Field)),string(CtrlVar.MustBe.(Field)))
    
    fprintf(" CtrlVar.%s has the value: \n \t  %s \n",Field,string(CtrlVar.(Field)))
    fprintf(" But this is not an allowed value! \n")
    fprintf(" The allowed values are: \n")
    fprintf('\t %s \n ',string(CtrlVar.MustBe.(Field)))
    error("Ua:CtrlVar","Invalid value")
end

end
