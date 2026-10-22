
function [InvertString,status]=SearchAndReplaceInverseFieldsInCtrlVar(InvertString)


% On return "InvertString" should have one of the following forms: 
%
%  "-?AGlen-"  "-?C-" , "-AGlen-?C-" 
%   where ? is either "" or "log"
%
%

status=true;

InvertString=string(InvertString) ;

InvertString=replace(InvertString,"Aglen","A") ; 

% first replace AGlen with A, although later I will return AGlen
% doing this replacement reduces the number of cases to look at
% and it is always possible to replace AGlen with A, whereas replacing A with
% AGlen is potentially ambiguous



switch lower(InvertString)
    
    case {"c","-c-"}
        
        InvertString="-C-";
        
    case {"a","-a-"}
        
        InvertString="-AGlen-";
        
    case {"b","-b-"}
        
        InvertString="-B-";
        
    case {"ac","ca","-a-c-","-c-a-","-ac-","-ca-","a-c","c-a"}
        
        InvertString="-AGlen-C-";
        
    case {"logc","-logc-"}
        
        InvertString="-logC-";
        
    case {"loga","-loga-"}
        
        InvertString="-logAGlen-";
        
    case {"-loga-logc-","-logalogc-","logalogc","loga-logc","-logc-loga-","-logcloga-","logcloga","logc-loga"}
        
        InvertString="-logAGlen-logC-";
        
    case {"-loga-c-","logac","loga-c","-logac-","-c-loga-","cloga","c-loga","-cloga-"}
        
        InvertString="-logAGlen-C-";
        
    case {"-a-logc-","alogc","a-logc","-alogc-","-logc-a-","logca","logc-a","-logca-"}
        
        InvertString="-AGlen-logC-";
        
    case {"-b-logc-","-logc-b-"}

        InvertString="-B-logC-";

    otherwise
        
        status=false;

end

%% This is a newer and simpler approach 
if ~status
    InvertStringOut="";

    if contains(InvertString,"logA")
        InvertStringOut=InvertStringOut+"-logAGlen-";
    end


    if contains(InvertString,"B")
        InvertStringOut=InvertStringOut+"-B-";
    end


    if contains(InvertString,"logC")
        InvertStringOut=InvertStringOut+"-logC-";
    end

    InvertStringOut=replace(InvertStringOut,"--","-");

    if InvertStringOut==""
        status=false;
    else
        InvertString=InvertStringOut ; 
        status=true;
    end


end




% does the string contain a but not aglen?  If so then something went wrong
if contains(lower(InvertString),"a") && ~contains(InvertString,"AGlen")
    status=false;
end




end