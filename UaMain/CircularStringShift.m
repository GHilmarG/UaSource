
function [NewString,FirstInstance]=CircularStringShift(String)


%% Get rid of eventual dashes at beginning and end
if String(1,1)=='-'
    String=String(1,2:end);
end


if String(1,end)=='-'
    String=String(1,1:end-1);
end

%%

NewString=String;
FirstInstance=String;
%% Now split and rejoin after a circular shift

temp=split(String,'-');

if numel(temp)>1
    NewString=join([temp(2:end) ; temp(1)],'-');
end

end