
function CtrlVar=ReplaceStructureFields(CtrlVar,CtrlVarOnInput)

Fields0=fieldnames(CtrlVarOnInput);
for I0 = 1:numel(Fields0)
    if isstruct(CtrlVarOnInput.(Fields0{I0}))
        Fields1=fieldnames(CtrlVarOnInput.(Fields0{I0})) ;
        for I1 = 1:numel(Fields1)
            if isstruct(CtrlVarOnInput.(Fields0{I0}).(Fields1{I1}))
                Fields2=fieldnames(CtrlVarOnInput.(Fields0{I0}).(Fields1{I1}));
                for I2 = 1:numel(Fields2)
                    if isstruct(CtrlVarOnInput.(Fields0{I0}).(Fields1{I1}).(Fields2{I2}))
                        error('only two sub structures allowed')
                    else
                        CtrlVar.(Fields0{I0}).(Fields1{I1}).(Fields2{I2})=CtrlVarOnInput.(Fields0{I0}).(Fields1{I1}).(Fields2{I2});
                    end
                end
            else
                CtrlVar.(Fields0{I0}).(Fields1{I1})=CtrlVarOnInput.(Fields0{I0}).(Fields1{I1});
            end
        end
    else
        CtrlVar.(Fields0{I0}) = CtrlVarOnInput.(Fields0{I0});
    end
end

end
