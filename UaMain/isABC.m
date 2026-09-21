




function [isA,isB,isC] = isABC(CtrlVar)




if contains(CtrlVar.Inverse.InvertFor,"-AGlen-")
    error("isABC:IncorrectInputs","Inversion for AGlen no longer supported. Invert for logAGlen instead.")
end


if contains(CtrlVar.Inverse.InvertFor,"-C-")
    error("isABC:IncorrectInputs","Inversion for C no longer supported. Invert for logC instead.")
end



if contains(CtrlVar.Inverse.InvertFor,"-logAGlen-")
    isA=true;
else
    isA=false;
end

if contains(CtrlVar.Inverse.InvertFor,'-B-')
    isB=true;
else 
    isB=false;
end

if contains(CtrlVar.Inverse.InvertFor,"-logC-")
    isC=true;
else
    isC=false; 
end



end