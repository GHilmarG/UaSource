
function [InvertString,status]=SearchAndReplaceInverseFieldsInCtrlVar(InvertString)

status=true;

switch lower(InvertString)
    
    case {'c','-c-'}
        InvertString='-C-';
    case {'aglen','-a-'}
        InvertString='-AGlen-';
    case {'b','-b-'}
        if contains(InvertString,'B')
            InvertString='-B-';
        elseif contains(InvertString,'b')
            error(' b inversion not possible, maybe try B inversion instead.') 
        end
    case {'aglenc','caglen','-aglen-c-','-c-aglen-','ac','ca','-a-c-','-c-a-'}
        InvertString='-AGlen-C-';
    case {'logc','-logc-'}
        InvertString='-logC-';
    case {'loga','logaglen','-loga-','-logaglen-'}
        InvertString='-logAGlen-';
    case {'logaglenlogc','logaglen-logc','-loga-logc-','-logc-loga-','logclogaglen','logc-logaglen','-logaglen-logc-','-logc-logaglen-'}
        InvertString='-logAGlen-logC-';
    otherwise
        
        status=false;
        
end

end