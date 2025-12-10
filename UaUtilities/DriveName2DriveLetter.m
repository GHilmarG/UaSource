function driveletter=DriveName2DriveLetter(drivename)


%
% driveletter=DriveName2DriveLetter(drivename)
%
% Finds the letter of the drive with the name `drivename'
%
% Example: driveletter=DriveName2DriveLetter('OS')  ; % on a windows system will most likely return  'C:\'
%
% Example
% cd(DriveName2DriveLetter('WorkSSD1'))  ; % will change to the top directory of the drive with the name `WorkSSD1'  (assuming there is a drive with such a name)
%

driveletter=[];
[driveletters,drivenames]=DriveLettersAndNames;



for n=1:numel(drivenames)
     if isequal(drivename,drivenames{n})
         driveletter=driveletters{n};
     end
end





end