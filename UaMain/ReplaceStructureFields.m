
function CtrlVar=ReplaceStructureFields(CtrlVar,CtrlVarOnInput)

%%
%
% Replaces structure fields in CtrlVar by those in CtrlVarOnInput
%
% Limitations:  Only up to three sub-structures are considered.
%
% Example:
%
%   CrlVar.a="a" ;
%   CtrlVar.b="b" ;
%   CtrlVar.c="c" ;
% 
%   CtrlVar.Field.a="a";
%   CtrlVar.Field.b="b";
%   CtrlVar.Field.c="c";
% 
%   CtrlVar.Field.Field.a="a";
%   CtrlVar.Field.Field.b="b";
% 
%   CtrlVar.Field.Field.Field.a="a";  
%   CtrlVar.Field.Field.Field.a="b";
%
%   CtrlVar.Field.Field.Field.Field.a="a";
% 
%   CtrlVarOnInput.a="A";
%   CtrlVarOnInput.c="C";
%   CtrlVarOnInput.Field.a="A";
%   CtrlVarOnInput.Field.Field.b="B";
%   CtrlVarOnInput.Field.Field.Field.a="A";  
%
% % CtrlVarOnInput.Field.Field.Field.Field.a="A";  % this should generate an error
%
%
%   CtrlVar=ReplaceStructureFields(CtrlVar,CtrlVarOnInput); 
%
%%


Fields0=fieldnames(CtrlVarOnInput);
for I0 = 1:numel(Fields0)
    if isstruct(CtrlVarOnInput.(Fields0{I0}))
        Fields1=fieldnames(CtrlVarOnInput.(Fields0{I0})) ;
        for I1 = 1:numel(Fields1)
            if isstruct(CtrlVarOnInput.(Fields0{I0}).(Fields1{I1}))
                Fields2=fieldnames(CtrlVarOnInput.(Fields0{I0}).(Fields1{I1}));
                for I2 = 1:numel(Fields2)
                    if isstruct(CtrlVarOnInput.(Fields0{I0}).(Fields1{I1}).(Fields2{I2}))
                        Fields3=fieldnames(CtrlVarOnInput.(Fields0{I0}).(Fields1{I1}).(Fields2{I2}));
                        for I3 = 1:numel(Fields3)
                            if isstruct(CtrlVarOnInput.(Fields0{I0}).(Fields1{I1}).(Fields2{I2}).(Fields3{I3}))
                                error('only three sub structures allowed')
                            else
                                CtrlVar.(Fields0{I0}).(Fields1{I1}).(Fields2{I2}).(Fields3{I3})=CtrlVarOnInput.(Fields0{I0}).(Fields1{I1}).(Fields2{I2}).(Fields3{I3});
                            end
                        end
                        
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
