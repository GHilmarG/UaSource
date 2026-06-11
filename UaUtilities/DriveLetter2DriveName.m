function name=DriveLetter2DriveName(DriveLetter)


% name=DriveLetter2DriveName(DriveLetter)
%
% finds name of a drive given the drive's letter
%
% example:  DriveName=DriveName('c')  ; % gives the name of the c drive
%

DriveLetter=lower(DriveLetter);
DriveLetter=DriveLetter(1:1);

[result,vol]=system(['vol ',DriveLetter,':']);

if result==0
    a=regexp(vol,'(?m)(?<drive>[A-Z]) is (?<label>\w+)$','names');
    if ~isempty(a)
    name=a.label;
    else 
        name=[];
    end
else
    name=[];
end


end